#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on 09.02.2015

 ----------------------------------------------------------------------------
* "THE BEER-WARE LICENSE" (Revision 42):
* <daniel.kratzert@uni-freiburg.de> wrote this file. As long as you retain this
* notice you can do whatever you want with this stuff. If we meet some day, and
* you think this stuff is worth it, you can buy me a beer in return.
* ----------------------------------------------------------------------------

@author: Daniel Kratzert
"""
from __future__ import print_function

from os.path import isfile

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from displaymol.sdm import SDM
from p4pfile.p4p_reader import P4PFile, read_file_to_list
from searcher.database_handler import Structure, get_cell_by_id, get_symmcards, get_atoms_table, \
    get_residuals, find_cell_by_volume, get_cells_as_list, get_all_structure_names, find_by_date, init_textsearch, \
    populate_fulltext_search_table, find_by_strings, find_by_it_number, find_by_elements, Base, get_lastrow_id
from searcher.filecrawler import fill_db_with_cif_data, put_files_in_db
from shelxfile.shelx import ShelXFile

DEBUG = True
import math
import os
import pathlib
import shutil
import sys
import tempfile
import time
from datetime import date

from PyQt5 import QtWidgets, QtCore, QtGui, uic

import misc.update_check
from apex import apeximporter
from displaymol import mol_file_writer, write_html
from lattice import lattice
from misc import update_check
from misc.version import VERSION
from pymatgen.core import mat_lattice
from searcher import constants, misc  # , database_handler
from searcher.constants import py36
from searcher.fileparser import Cif
from searcher.misc import is_valid_cell, elements
from contextlib import contextmanager

if py36:
    """Only import this if Python 3.6 is used."""
    try:
        from PyQt5.QtWebEngineWidgets import QWebEngineView
    except Exception as e:
        print(e)

__metaclass__ = type  # use new-style classes

"""
TODO:
- Improve text search (in cif file). Figure out which tokenchars configuration works best.
- refractor put_cifs_in_db()
- Move indexer to extra thread
- sort results by G6 distance
- colour table according to structure quality.
  
Search for:
- draw structure (with JSME? Acros? Kekule?, https://github.com/ggasoftware/ketcher)
- compare  molecules https://groups.google.com/forum/#!msg/networkx-discuss/gC_-Wc0bRWw/ISRZYFsPCQAJ
  - search algorithms
  http://chemmine.ucr.edu/help/#similarity, https://en.wikipedia.org/wiki/Jaccard_index
"""
# This is to make sure that strf finds the application path even when it is
# executed from another path e.g. when opened via "open file" in windows:
if getattr(sys, 'frozen', False):
    # If the application is run as a bundle, the pyInstaller bootloader
    # extends the sys module by a flag frozen=True and sets the app
    # path into variable _MEIPASS'.
    application_path = sys._MEIPASS
else:
    application_path = os.path.dirname(os.path.abspath(__file__))


class StartStructureDB(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ui = Ui_stdbMainwindow()
        self.ui.setupUi(self)
        self.statusBar().showMessage('StructureFinder version {}'.format(VERSION))
        self.ui.cifList_treeWidget.show()
        self.ui.cifList_treeWidget.hideColumn(3)
        self.dbfdesc = None
        self.dbfilename = None
        self.tmpfile = False  # indicates wether a tmpfile or any other db file is used
        # self.ui.centralwidget.setMinimumSize(1000, 500)
        self.abort_import_button = QtWidgets.QPushButton("Abort", parent=self.ui.statusbar)
        self.progress = QtWidgets.QProgressBar(self)
        self.progress.setFormat('')
        self.ui.statusbar.addWidget(self.progress)
        self.ui.statusbar.addWidget(self.abort_import_button)
        self.Session = None
        self.apx = None
        self.structureId = '0'
        self.passwd = ''
        self.show()
        self.setAcceptDrops(True)
        self.full_list = True  # indicator if the full structures list is shown
        self.decide_import = True
        self.connect_signals_and_slots()
        # Set both to today() to distinquish between a modified and unmodified date field.
        self.ui.dateEdit1.setDate(QtCore.QDate(date.today()))
        self.ui.dateEdit2.setDate(QtCore.QDate(date.today()))
        if py36:
            try:
                molf = pathlib.Path(os.path.join(application_path, "./displaymol/jsmol.htm"))
                molf.write_text(data=' ', encoding="utf-8", errors='ignore')
                self.init_webview()
            except Exception as e:
                # Graphics driver not compatible
                print(e)
        else:
            self.ui.tabWidget.removeTab(2)
            self.ui.txtSearchEdit.hide()
            self.ui.txtSearchLabel.hide()
            self.ui.openglview.hide()
        self.ui.tabWidget.setCurrentIndex(0)
        self.setWindowIcon(QtGui.QIcon(os.path.join(application_path, './icons/strf.png')))
        self.uipass = Ui_PasswdDialog()
        # self.ui.cifList_treeWidget.sortByColumn(0, 0)
        # Actions for certain gui elements:
        self.ui.cellField.addAction(self.ui.actionCopy_Unit_Cell)
        self.ui.cifList_treeWidget.addAction(self.ui.actionGo_to_All_CIF_Tab)
        if len(sys.argv) > 1:
            self.dbfilename = sys.argv[1]
            if isfile(self.dbfilename):
                self.engine = create_engine('sqlite:///' + self.dbfilename)
                self.Session = sessionmaker(bind=self.engine)
                self.show_full_list()
        if update_check.is_update_needed(VERSION=VERSION):
            self.statusBar().showMessage('A new Version of StructureFinder is available at '
                                         'https://www.xs3.uni-freiburg.de/research/structurefinder')
        # select the first item in the list
        item = self.ui.cifList_treeWidget.topLevelItem(0)
        self.ui.cifList_treeWidget.setCurrentItem(item)

    @contextmanager
    def session_scope(self):
        """Provide a transactional scope around a series of operations."""
        session = self.Session()
        # session.autoflush = False
        try:
            yield session
            session.commit()
        except:
            session.rollback()
            raise
        finally:
            session.close()

    def connect_signals_and_slots(self):
        """
        Connects the signals and slot.
        The actionExit signal is connected in the ui file.
        """
        # Buttons:
        self.ui.importDatabaseButton.clicked.connect(self.open_sql_database)
        self.ui.saveDatabaseButton.clicked.connect(self.save_database)
        self.ui.importDirButton.clicked.connect(self.import_cif_dirs)
        self.ui.openApexDBButton.clicked.connect(self.import_apex_db)
        self.ui.closeDatabaseButton.clicked.connect(self.close_db)
        self.abort_import_button.clicked.connect(self.abort_import)
        self.ui.moreResultsCheckBox.stateChanged.connect(self.cell_state_changed)
        self.ui.sublattCheckbox.stateChanged.connect(self.cell_state_changed)
        self.ui.ad_SearchPushButton.clicked.connect(self.advanced_search)
        self.ui.ad_ClearSearchButton.clicked.connect(self.show_full_list)
        # Actions:
        self.ui.actionClose_Database.triggered.connect(self.close_db)
        self.ui.actionImport_directory.triggered.connect(self.import_cif_dirs)
        self.ui.actionImport_file.triggered.connect(self.open_sql_database)
        self.ui.actionSave_Database.triggered.connect(self.save_database)
        self.ui.actionCopy_Unit_Cell.triggered.connect(self.copyUnitCell)
        self.ui.actionGo_to_All_CIF_Tab.triggered.connect(self.on_click_item)
        # Other fields:
        if py36:
            self.ui.txtSearchEdit.textChanged.connect(self.search_text)
        else:
            self.ui.txtSearchEdit.setText("For full test search, use a modern Operating system.")
        self.ui.searchCellLineEDit.textChanged.connect(self.search_cell)
        self.ui.p4pCellButton.clicked.connect(self.get_name_from_p4p)
        self.ui.cifList_treeWidget.selectionModel().currentChanged.connect(self.get_properties)
        self.ui.cifList_treeWidget.itemDoubleClicked.connect(self.on_click_item)
        self.ui.ad_elementsIncLineEdit.textChanged.connect(self.elements_fields_check)
        self.ui.ad_elementsExclLineEdit.textChanged.connect(self.elements_fields_check)
        self.ui.add_res.clicked.connect(self.res_checkbox_clicked)
        self.ui.add_cif.clicked.connect(self.cif_checkbox_clicked)
        self.ui.growCheckBox.toggled.connect(self.redraw_molecule)

    def res_checkbox_clicked(self, click):
        if not any([self.ui.add_res.isChecked(), self.ui.add_cif.isChecked()]):
            self.ui.add_cif.setChecked(True)

    def cif_checkbox_clicked(self, click):
        if not any([self.ui.add_res.isChecked(), self.ui.add_cif.isChecked()]):
            self.ui.add_cif.setChecked(True)

    def on_click_item(self, item):
        self.ui.tabWidget.setCurrentIndex(1)

    def dragEnterEvent(self, e):
        if e.mimeData().hasText():
            e.accept()
        else:
            e.ignore()

    def dropEvent(self, e):
        """
        Handles drop events. 
        """
        from urllib.parse import urlparse
        p = urlparse(e.mimeData().text())
        if sys.platform.startswith('win'):
            final_path = p.path[1:]  # remove strange / at start
        else:
            final_path = p.path
        _, ending = os.path.splitext(final_path)
        # print(final_path, ending)
        if ending == '.p4p':
            self.search_for_p4pcell(final_path)
        if ending == '.res' or ending == '.ins':
            self.search_for_res_cell(final_path)
        if ending == '.cif':
            self.search_for_cif_cell(final_path)

    @staticmethod
    def validate_sumform(inelem: list):
        """
        Checks if the elements in inelem are valid Chemical elements
        """
        ok = True
        for el in inelem:
            if el not in elements:
                ok = False
        return ok

    @QtCore.pyqtSlot('QString', name="elements_fields_check")
    def elements_fields_check(self):
        """
        """
        elem1 = self.ui.ad_elementsIncLineEdit.text().split()
        elem2 = self.ui.ad_elementsExclLineEdit.text().split()
        if (not self.elements_doubled_check(elem1, elem2)) or (not self.elements_doubled_check(elem2, elem1)):
            self.elements_invalid()
        else:
            self.elements_regular()

    def elements_doubled_check(self, elem1, elem2):
        """
        Validates if elements of elem1 are not in elem2 and elem1 has purely valid sum formula. 
        """
        ok = True
        for el in elem1:
            if el in elem2:
                ok = False
        if not self.validate_sumform(elem1):
            ok = False
        return ok

    def elements_invalid(self):
        # Elements not valid:
        self.ui.ad_elementsIncLineEdit.setStyleSheet("color: rgb(255, 0, 0); font: bold 12px;")
        self.ui.ad_SearchPushButton.setDisabled(True)
        self.ui.ad_elementsExclLineEdit.setStyleSheet("color: rgb(255, 0, 0); font: bold 12px;")
        self.ui.ad_SearchPushButton.setDisabled(True)

    def elements_regular(self):
        # Elements valid:
        self.ui.ad_elementsIncLineEdit.setStyleSheet("color: rgb(0, 0, 0);")
        self.ui.ad_SearchPushButton.setEnabled(True)
        self.ui.ad_elementsExclLineEdit.setStyleSheet("color: rgb(0, 0, 0);")
        self.ui.ad_SearchPushButton.setEnabled(True)

    def copyUnitCell(self):
        if self.structureId:
            try:
                with self.session_scope() as session:
                    cell = "{:>6.3f} {:>6.3f} {:>6.3f} {:>6.3f} {:>6.3f} {:>6.3f}" \
                        .format(*get_cell_by_id(session, self.structureId))
            except Exception:
                return False
            clipboard = QtWidgets.QApplication.clipboard()
            clipboard.setText(cell)
            self.ui.statusbar.showMessage('Copied unit cell {} to clip board.'
                                          .format(cell))
        return True

    @QtCore.pyqtSlot(name="advanced_search")
    def advanced_search(self):
        """
        Combines all the search fields. Collects all includes, all excludes ad calculates
        the difference.
        """
        excl = []
        incl = []
        date_results = []
        results = []
        it_results = []
        cell = is_valid_cell(self.ui.ad_unitCellLineEdit.text())
        date1 = self.ui.dateEdit1.text()
        date2 = self.ui.dateEdit2.text()
        elincl = self.ui.ad_elementsIncLineEdit.text().strip(' ')
        elexcl = self.ui.ad_elementsExclLineEdit.text().strip(' ')
        txt = self.ui.ad_textsearch.text().strip(' ')
        txt_ex = self.ui.ad_textsearch_excl.text().strip(' ')
        spgr = self.ui.SpGrcomboBox.currentText()
        try:
            spgr = int(spgr.split()[0])
        except:
            spgr = 0
        if spgr:
            with self.session_scope() as session:
                it_results = find_by_it_number(session, spgr)
        if date1 != date2:
            date_results = self.find_dates(date1, date2)
        if cell:
            cellres = self.search_cell_idlist(cell)
            incl.append(cellres)
        if elincl:
            incl.append(self.search_elements(elincl))
        if txt:
            if len(txt) >= 2 and "*" not in txt:
                txt = '*' + txt + '*'
            idlist = [i[0] for i in find_by_strings(self.engine, txt)]
            if idlist:
                incl.append(idlist)
        if elexcl:
            excl.append(self.search_elements(elexcl, anyresult=True))
        if txt_ex:
            if len(txt_ex) >= 2 and "*" not in txt_ex:
                txt_ex = '*' + txt_ex + '*'
            idlist = [i[0] for i in find_by_strings(self.engine, txt_ex)]
            if idlist:
                excl.append(idlist)
        if incl and incl[0]:
            results = set(incl[0]).intersection(*incl)
            if date_results:
                results = set(date_results).intersection(results)
            if it_results:
                results = set(it_results).intersection(results)
        elif date_results and not it_results:
            results = date_results
        elif not date_results and it_results:
            results = it_results
        elif it_results and date_results:
            results = set(it_results).intersection(date_results)
        if not results:
            self.statusBar().showMessage('Found 0 structures.')
            return
        if excl:
            try:
                self.display_structures_by_idlist(list(results - set(misc.flatten(excl))))
            except TypeError:
                return
        else:
            self.display_structures_by_idlist(list(results))

    def display_structures_by_idlist(self, idlist: list or set) -> None:
        """
        Displays the structures with id in results list
        """
        if not idlist:
            self.statusBar().showMessage('Found {} structures.'.format(0))
            return
        with self.session_scope() as session:
            searchresult = get_all_structure_names(session, idlist)
        self.statusBar().showMessage('Found {} structures.'.format(len(idlist)))
        self.ui.cifList_treeWidget.clear()
        self.full_list = False
        for row in searchresult:
            self.add_table_row(filename=row.filename, path=row.path, structure_id=row.Id, data=row.dataname)
        self.set_columnsize()
        # self.ui.cifList_treeWidget.resizeColumnToContents(0)
        if idlist:
            self.ui.tabWidget.setCurrentIndex(0)

    def passwd_handler(self):
        """
        Manages the QDialog buttons on the password dialog.
        """
        d = QtWidgets.QDialog()
        self.passwd = self.uipass.setupUi(d)
        ip_dialog = d.exec()
        if ip_dialog == 1:  # Accepted
            self.import_apex_db(user=self.uipass.userNameLineEdit.text(),
                                password=self.uipass.PasswordLineEdit.text(),
                                host=self.uipass.IPlineEdit.text())
        elif ip_dialog == 0:  # reject
            return None
        else:
            return None

    def init_webview(self):
        """
        Initializes a QWebengine to view the molecule.
        """
        self.view = QWebEngineView()
        # QtWebEngine.initialize()
        self.view.load(
                QtCore.QUrl.fromLocalFile(os.path.abspath(os.path.join(application_path, "./displaymol/jsmol.htm"))))
        # self.view.setMaximumWidth(260)
        # self.view.setMaximumHeight(290)
        self.ui.ogllayout.addWidget(self.view)
        # self.view.show()

    @QtCore.pyqtSlot(name="cell_state_changed")
    def cell_state_changed(self):
        """
        Searches a cell but with diffeent loose or strict option.
        """
        self.search_cell(self.ui.searchCellLineEDit.text())

    def import_cif_dirs(self):
        # worker = RunIndexerThread(self)
        # worker.start()
        self.tmpfile = True
        self.statusBar().showMessage('')
        self.close_db()
        self.start_db()
        self.progressbar(1, 0, 20)
        self.abort_import_button.show()
        fname = QtWidgets.QFileDialog.getExistingDirectory(self, 'Open Directory', '')
        if not fname:
            self.progress.hide()
            self.abort_import_button.hide()
        with self.session_scope() as session:
            put_files_in_db(self, searchpath=fname, fillres=self.ui.add_res.isChecked(),
                            fillcif=self.ui.add_cif.isChecked(), session=session)
        self.progress.hide()
        init_textsearch(self.engine)
        populate_fulltext_search_table(self.engine)
        self.ui.cifList_treeWidget.show()
        self.set_columnsize()
        # self.ui.cifList_treeWidget.resizeColumnToContents(0)
        # self.ui.cifList_treeWidget.resizeColumnToContents(1)
        # self.ui.cifList_treeWidget.sortByColumn(0, 0)
        self.abort_import_button.hide()

    def progressbar(self, curr: int, min: int, max: int) -> None:
        """
        Displays a progress bar in the status bar.
        """
        self.progress.setValue(curr)
        self.progress.setMaximum(max)
        self.progress.setMinimum(min)
        self.progress.show()
        if curr == max:
            self.progress.hide()

    @QtCore.pyqtSlot(name="close_db")
    def close_db(self, copy_on_close: str = None) -> bool:
        """
        Closed the current database and erases the list.
        :param copy_on_close: Path to where the file should be copied after close()
        """
        self.ui.searchCellLineEDit.clear()
        self.ui.txtSearchEdit.clear()
        self.ui.cifList_treeWidget.clear()
        if py36:
            molf = pathlib.Path(os.path.join(application_path, "./displaymol/jsmol.htm"))
            molf.write_text(data=' ', encoding="utf-8", errors='ignore')
            self.view.reload()
        if copy_on_close:
            if shutil._samefile(self.dbfilename, copy_on_close):
                self.statusBar().showMessage("You can not save to the currently opened file!", msecs=5000)
                return False
            else:
                shutil.copy(self.dbfilename, copy_on_close)
            if self.tmpfile:
                try:
                    os.remove(self.dbfilename)
                    self.dbfilename = None
                except Exception:
                    return False
        return True

    @QtCore.pyqtSlot(name="abort_import")
    def abort_import(self):
        """
        This slot means, import was aborted.
        """
        self.decide_import = False

    def start_db(self):
        """
        Initializes the database.
        """
        self.dbfdesc, self.dbfilename = tempfile.mkstemp()
        self.engine = create_engine('sqlite:///' + self.dbfilename)
        # engine.echo = True
        self.Session = sessionmaker(bind=self.engine)
        Base.metadata.create_all(self.engine)  # creates the table

    @QtCore.pyqtSlot('QModelIndex', name="get_properties")
    def get_properties(self, item):
        """
        This slot shows the properties of a cif file in the properties widget
        """
        structure_id = item.sibling(item.row(), 3).data()
        with self.session_scope() as session:
            dic = get_residuals(session, structure_id)
        if not dic:
            return False
        self.structureId = structure_id
        self.display_properties(structure_id, dic)
        self.structureId = structure_id
        return True

    def save_database(self) -> bool:
        """
        Saves the database to a certain file. Therefore I have to close the database.
        """
        status = False
        save_name, tst = QtWidgets.QFileDialog.getSaveFileName(self, caption='Save File', directory='./',
                                                               filter="*.sqlite")
        if save_name:
            if shutil._samefile(self.dbfilename, save_name):
                self.statusBar().showMessage("You can not save to the currently opened file!", msecs=5000)
                return False
        if save_name:
            status = self.close_db(save_name)
        if status:
            self.statusBar().showMessage("Database saved.", msecs=5000)

    def eventFilter(self, object, event):
        """Event filter for mouse clicks."""
        if event.type() == QtCore.QEvent.MouseButtonDblClick:
            self.copyUnitCell()
        elif event.type() == QtCore.QEvent.MouseButtonPress:
            if event.buttons() == QtCore.Qt.RightButton:
                # print("rightbutton")
                return True
        return False

    def redraw_molecule(self):
        with self.session_scope() as session:
            cell = get_cell_by_id(session, self.structureId)
        if not cell:
            return False
        try:
            self.display_molecule(cell, str(self.structureId))
        except Exception as e:
            print(e, "unable to display molecule")
            if DEBUG:
                raise

    def display_properties(self, structure_id: str, cif_dic: dict) -> bool:
        """
        Displays the residuals from the cif file
        Measured Refl.
        Independent Refl.
        Reflections Used

        _refine_ls_number_reflns -> unique reflect. (Independent reflections)
        _reflns_number_gt        -> unique über 2sigma (Independent reflections >2sigma)
        """
        self.clear_fields()
        with self.session_scope() as session:
            cell = get_cell_by_id(session, structure_id)
        a, b, c, alpha, beta, gamma, volume = cell
        if not all((a, b, c, alpha, beta, gamma, volume)):
            self.ui.cellField.setText('            ')
            return False
        if not cell:
            return False
        try:
            self.display_molecule(cell[:6], structure_id)
        except Exception as e:
            print(e, "unable to display molecule")
            if DEBUG:
                raise
        self.ui.cifList_treeWidget.setFocus()
        if not cif_dic:
            return False
        # self.ui.cellField.setMinimumWidth(180)
        self.ui.cellField.setText(constants.celltxt.format(a, alpha, b, beta, c, gamma, volume, ''))
        self.ui.cellField.installEventFilter(self)
        self.ui.cellField.setToolTip("Double click on 'Unit Cell' to copy to clipboard.")
        try:
            self.ui.wR2LineEdit.setText("{:>5.4f}".format(cif_dic['_refine_ls_wR_factor_ref']))
        except (ValueError, TypeError):
            pass
        try:  # R1:
            if cif_dic['_refine_ls_R_factor_gt']:
                self.ui.r1LineEdit.setText("{:>5.4f}".format(cif_dic['_refine_ls_R_factor_gt']))
            else:
                self.ui.r1LineEdit.setText("{:>5.4f}".format(cif_dic['_refine_ls_R_factor_all']))
        except (ValueError, TypeError):
            pass
        self.ui.zLineEdit.setText("{}".format(cif_dic['_cell_formula_units_Z']))
        try:
            sumform = misc.format_sum_formula(cif_dic['_chemical_formula_sum'])
        except KeyError:
            sumform = ''
        self.ui.formLabel.setText("{}".format(sumform))
        self.ui.reflTotalLineEdit.setText("{}".format(cif_dic['_diffrn_reflns_number']))
        self.ui.uniqReflLineEdit.setText("{}".format(cif_dic['_refine_ls_number_reflns']))
        self.ui.refl2sigmaLineEdit.setText("{}".format(cif_dic['_reflns_number_gt']))
        self.ui.goofLineEdit.setText("{}".format(cif_dic['_refine_ls_goodness_of_fit_ref']))
        it_num = cif_dic['_space_group_IT_number']
        if it_num:
            it_num = "({})".format(it_num)
        self.ui.SpaceGroupLineEdit.setText("{} {}".format(cif_dic['_space_group_name_H_M_alt'], it_num))
        self.ui.temperatureLineEdit.setText("{}".format(cif_dic['_diffrn_ambient_temperature']))
        self.ui.maxShiftLineEdit.setText("{}".format(cif_dic['_refine_ls_shift_su_max']))
        peak = cif_dic['_refine_diff_density_max']
        if peak:
            self.ui.peakLineEdit.setText("{} / {}".format(peak, cif_dic['_refine_diff_density_min']))
        self.ui.rintLineEdit.setText("{}".format(cif_dic['_diffrn_reflns_av_R_equivalents']))
        self.ui.rsigmaLineEdit.setText("{}".format(cif_dic['_diffrn_reflns_av_unetI_netI']))
        self.ui.cCDCNumberLineEdit.setText("{}".format(cif_dic['_database_code_depnum_ccdc_archive']))
        try:
            dat_param = cif_dic['_refine_ls_number_reflns'] / cif_dic['_refine_ls_number_parameters']
        except (ValueError, ZeroDivisionError, TypeError):
            dat_param = 0.0
        self.ui.dataReflnsLineEdit.setText("{:<5.1f}".format(dat_param))
        self.ui.numParametersLineEdit.setText("{}".format(cif_dic['_refine_ls_number_parameters']))
        wavelen = cif_dic['_diffrn_radiation_wavelength']
        thetamax = cif_dic['_diffrn_reflns_theta_max']
        # d = lambda/2sin(theta):
        try:
            d = wavelen / (2 * math.sin(math.radians(thetamax)))
        except(ZeroDivisionError, TypeError):
            d = 0.0
        self.ui.numRestraintsLineEdit.setText("{}".format(cif_dic['_refine_ls_number_restraints']))
        self.ui.thetaMaxLineEdit.setText("{}".format(thetamax))
        self.ui.thetaFullLineEdit.setText("{}".format(cif_dic['_diffrn_reflns_theta_full']))
        self.ui.dLineEdit.setText("{:5.3f}".format(d))
        self.ui.lastModifiedLineEdit.setText(str(cif_dic['modification_time']))
        try:
            compl = cif_dic['_diffrn_measured_fraction_theta_max'] * 100
            if not compl:
                compl = 0.0
        except TypeError:
            compl = 0.0
        try:
            self.ui.completeLineEdit.setText("{:<5.1f}".format(compl))
        except ValueError:
            pass
        self.ui.wavelengthLineEdit.setText("{}".format(wavelen))
        #####  Maybe put this into a new method and only let it run when all values tab opens:
        self.ui.allCifTreeWidget.clear()
        # This makes selection slow and is not really needed:
        # atoms_item = QtWidgets.QTreeWidgetItem()
        # self.ui.allCifTreeWidget.addTopLevelItem(atoms_item)
        # atoms_item.setText(0, 'Atoms')
        # try:
        #    atoms = self.structures.get_atoms_table(structure_id, cartesian=False)
        #    for at in atoms:
        #        data_cif_tree_item = QtWidgets.QTreeWidgetItem(atoms_item)
        #        self.ui.allCifTreeWidget.addTopLevelItem(atoms_item)
        #        data_cif_tree_item.setText(1, '{:<8.8s}\t {:<4s}\t {:>8.5f}\t {:>8.5f}\t {:>8.5f}'.format(*at))
        # except TypeError:
        #    pass
        #######################################
        for key, value in cif_dic.items():
            if key == "_shelx_res_file":
                self.ui.SHELXplainTextEdit.setPlainText(cif_dic['_shelx_res_file'])
                continue
            cif_tree_item = QtWidgets.QTreeWidgetItem()
            self.ui.allCifTreeWidget.addTopLevelItem(cif_tree_item)
            cif_tree_item.setText(0, str(key))
            cif_tree_item.setText(1, str(value))
        if not cif_dic['_shelx_res_file']:
            self.ui.SHELXplainTextEdit.setPlainText("No SHELXL res file in cif found.")
        # self.ui.cifList_treeWidget.sortByColumn(0, 0)
        self.ui.allCifTreeWidget.resizeColumnToContents(0)
        self.ui.allCifTreeWidget.resizeColumnToContents(1)
        return True

    def display_molecule(self, cell: list, structure_id: str) -> None:
        """
        Creates a html file from a mol file to display the molecule in jsmol-lite
        """
        with self.session_scope() as session:
            symmcards = get_symmcards(session, structure_id)
        grow = self.ui.growCheckBox.isChecked()
        if len(symmcards) == 1:
            grow = False
        blist = None
        if grow:
            with self.session_scope() as session:
                atoms = get_atoms_table(session, structure_id, cell, cartesian=False)
            if atoms:
                sdm = SDM(atoms, symmcards, cell)
                needsymm = sdm.calc_sdm()
                atoms = sdm.packer(sdm, needsymm)
        else:
            with self.session_scope() as session:
                atoms = get_atoms_table(session, structure_id, cell, cartesian=True)
        try:
            tst = mol_file_writer.MolFile(atoms, blist)
            mol = tst.make_mol()
        except (TypeError, KeyError) as e:
            print("Error in structure", structure_id, "while writing mol file.\n", e)
            mol = ' '
            if DEBUG:
                raise
        # print(self.ui.openglview.width()-30, self.ui.openglview.height()-50)
        content = write_html.write(mol, self.ui.openglview.width() - 30, self.ui.openglview.height() - 50)
        p2 = pathlib.Path(os.path.join(application_path, "./displaymol/jsmol.htm"))
        p2.write_text(data=content, encoding="utf-8", errors='ignore')
        self.view.reload()

    @QtCore.pyqtSlot('QString')
    def find_dates(self, date1: str, date2: str) -> list:
        """
        Returns a list if id between date1 and date2
        """
        if not date1:
            date1 = '0000-01-01'
        if not date2:
            date2 = 'NOW'
        with self.session_scope() as session:
            result = find_by_date(session, date1, date2)
        return result

    @QtCore.pyqtSlot('QString')
    def search_text(self, search_string: str) -> bool:
        """
        searches db for given text
        """
        self.ui.searchCellLineEDit.clear()
        self.ui.cifList_treeWidget.clear()
        result_list = []
        if len(search_string) == 0:
            self.show_full_list()
            return False
        if len(search_string) >= 2 and "*" not in search_string:
            search_string = "{}{}{}".format('*', search_string, '*')
        try:
            result_list = find_by_strings(self.engine, search_string)
        except AttributeError as e:
            print(e)
        try:
            self.statusBar().showMessage("Found {} entries.".format(len(result_list)))
            for row in result_list:
                self.add_table_row(filename=row[1], path=row[3], structure_id=row[0], data=row[2])
            self.set_columnsize()
            # self.ui.cifList_treeWidget.resizeColumnToContents(0)
        except Exception:
            self.statusBar().showMessage("Nothing found.")

    def search_cell_idlist(self, cell: list) -> list:
        """
        Searches for a unit cell and resturns a list of found database ids.
        This method does not validate the cell. This has to be done before!
        """
        if self.ui.moreResultsCheckBox.isChecked() or self.ui.ad_moreResultscheckBox.isChecked():
            # more results:
            vol_threshold = 0.09
            ltol = 0.2
            atol = 2
        else:
            # regular:
            vol_threshold = 0.03
            ltol = 0.06
            atol = 1
        try:
            volume = lattice.vol_unitcell(*cell)
            with self.session_scope() as session:
                idlist = find_cell_by_volume(session, volume, vol_threshold)
            if self.ui.sublattCheckbox.isChecked() or self.ui.ad_superlatticeCheckBox.isChecked():
                # sub- and superlattices:
                for v in [volume * x for x in [2.0, 3.0, 4.0, 6.0, 8.0, 10.0]]:
                    # First a list of structures where the volume is similar:
                    with self.session_scope() as session:
                        idlist.extend(find_cell_by_volume(session, v, vol_threshold))
                idlist = list(set(idlist))
                idlist.sort()
        except (ValueError, AttributeError):
            if not self.full_list:
                self.ui.cifList_treeWidget.clear()
                self.statusBar().showMessage('Found 0 cells.')
            return []
        # Real lattice comparing in G6:
        idlist2 = []
        if idlist:
            lattice1 = mat_lattice.Lattice.from_parameters_niggli_reduced(*cell)
            self.statusBar().clearMessage()
            cells = []
            # SQLite can only handle 999 variables at once:
            # for cids in chunks(idlist, 500):
            with self.session_scope() as session:
                cells.extend(get_cells_as_list(session, idlist))
            for num, cell_id in enumerate(idlist):
                self.progressbar(num, 0, len(idlist) - 1)
                try:
                    lattice2 = mat_lattice.Lattice.from_parameters(
                            float(cells[num][0]),
                            float(cells[num][1]),
                            float(cells[num][2]),
                            float(cells[num][3]),
                            float(cells[num][4]),
                            float(cells[num][5]))
                except ValueError:
                    continue
                mapping = lattice1.find_mapping(lattice2, ltol, atol, skip_rotation_matrix=True)
                if mapping:
                    # pprint.pprint(map[3])
                    idlist2.append(cell_id)
        # print("After match: ", len(idlist2))
        return idlist2

    @QtCore.pyqtSlot('QString', name='search_cell')
    def search_cell(self, search_string: str) -> bool:
        """
        searches db for given cell via the cell volume
        """
        cell = is_valid_cell(search_string)
        self.ui.txtSearchEdit.clear()
        if not cell:
            if self.ui.searchCellLineEDit.text():
                self.statusBar().showMessage('Not a valid unit cell!', msecs=3000)
                return False
            else:
                self.full_list = True  # Set status where full list is displayed
                self.show_full_list()
            if self.full_list:
                return False
            return False
        idlist = self.search_cell_idlist(cell)
        if not idlist:
            self.ui.cifList_treeWidget.clear()
            self.statusBar().showMessage('Found 0 cells.', msecs=0)
            return False
        self.ui.cifList_treeWidget.clear()
        with self.session_scope() as session:
            for row in get_all_structure_names(session, idlist):
                self.add_table_row(filename=row.filename, path=row.path, structure_id=row.Id, data=row.dataname)
        self.statusBar().showMessage('Found {} cells.'.format(len(idlist)))

        self.full_list = False
        self.set_columnsize()
        # self.ui.cifList_treeWidget.sortByColumn(0, 0)
        # self.ui.cifList_treeWidget.resizeColumnToContents(0)
        return True

    def search_elements(self, elements: str, anyresult: bool = False) -> list:
        """
        list(set(l).intersection(l2))
        """
        self.statusBar().showMessage('')
        res = []
        try:
            formula = misc.get_list_of_elements(elements)
        except KeyError:
            self.statusBar().showMessage('Error: Wrong list of Elements!', msecs=5000)
            return []
        try:
            with self.session_scope() as session:
                res = find_by_elements(session, formula, anyresult=anyresult)
        except AttributeError:
            pass
        return list(res)

    def add_table_row(self, filename: str, path: str, data: bytes, structure_id: str) -> None:
        """
        Adds a line to the search results table.
        """
        if isinstance(filename, bytes):
            filename = filename.decode("utf-8", "surrogateescape")
        if isinstance(path, bytes):
            path = path.decode("utf-8", "surrogateescape")
        if isinstance(data, bytes):
            data = data.decode("utf-8", "surrogateescape")
        tree_item = QtWidgets.QTreeWidgetItem()
        tree_item.setText(0, filename)  # name
        tree_item.setText(1, data)  # data
        tree_item.setText(2, path)  # path
        tree_item.setData(3, 0, structure_id)  # id
        self.ui.cifList_treeWidget.addTopLevelItem(tree_item)

    def open_sql_database(self) -> bool:
        """
        Import a new database.
        """
        self.tmpfile = False
        self.close_db()
        fname = QtWidgets.QFileDialog.getOpenFileName(self, caption='Open File', directory='./', filter="*.sqlite")
        if not fname[0]:
            return False
        print("Opened {}.".format(fname[0]))
        self.dbfilename = fname[0]
        try:
            self.engine.dispose()
        except AttributeError:
            pass
        self.engine = create_engine('sqlite:///' + self.dbfilename)
        if self.Session:
            self.Session.configure(bind=self.engine)
        else:
            self.Session = sessionmaker(bind=self.engine)
        self.show_full_list()
        if not self.Session:
            return False
        return True

    def open_apex_db(self, user: str, password: str, host: str) -> bool:
        """
        Opens the APEX db to be displayed in the treeview.
        """
        self.apx = apeximporter.ApexDB()
        connok = False
        try:
            connok = self.apx.initialize_db(user, password, host)
        except Exception:
            self.passwd_handler()
        return connok

    def get_name_from_p4p(self):
        """
        Reads a p4p file to get the included unit cell for a cell search.
        """
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, caption='Open p4p File', directory='./',
                                                         filter="*.p4p *.cif *.res *.ins")
        _, ending = os.path.splitext(fname)
        if ending == '.p4p':
            self.search_for_p4pcell(fname)
        if ending in ['.res', '.ins']:
            self.search_for_res_cell(fname)
        if ending == '.cif':
            self.search_for_cif_cell(fname)

    def search_for_p4pcell(self, fname):
        if fname:
            p4plist = read_file_to_list(fname)
            p4p = P4PFile(p4plist)
        else:
            return
        if p4p:
            if p4p.cell:
                try:
                    self.ui.searchCellLineEDit.setText('{:<6.3f} {:<6.3f} {:<6.3f} '
                                                       '{:<6.3f} {:<6.3f} {:<6.3f}'.format(*p4p.cell))
                except TypeError:
                    pass
            else:
                self.moving_message('Could not read P4P file!')
        else:
            self.moving_message('Could not read P4P file!')

    def search_for_res_cell(self, fname):
        if fname:
            shx = ShelXFile(fname)
        else:
            return
        if shx:
            if shx.cell:
                try:
                    self.ui.searchCellLineEDit.setText('{:<6.3f} {:<6.3f} {:<6.3f} '
                                                       '{:<6.3f} {:<6.3f} {:<6.3f}'.format(*shx.cell.cell_list))
                except TypeError:
                    pass
            else:
                self.moving_message('Could not read res file!')
        else:
            self.moving_message('Could not read res file!')

    def search_for_cif_cell(self, fname):
        if fname:
            cif = Cif()
            try:
                cif.parsefile(pathlib.Path(fname).read_text(encoding='utf-8',
                                                            errors='ignore').splitlines(keepends=True))
            except FileNotFoundError:
                self.moving_message('File not found.')
        else:
            return
        if cif:
            if cif.cell:
                try:
                    self.ui.searchCellLineEDit.setText('{:<6.3f} {:<6.3f} {:<6.3f} '
                                                       '{:<6.3f} {:<6.3f} {:<6.3f}'.format(*cif.cell[:6]))
                except TypeError:
                    pass
            else:
                self.moving_message('Could not read cif file!')
        else:
            self.moving_message('Could not read cif file!')

    def moving_message(self, message="", times=20):
        for s in range(times):
            time.sleep(0.05)
            self.statusBar().showMessage("{}{}".format(' ' * s, message))

    def import_apex_db(self, user: str = '', password: str = '', host: str = '') -> None:
        """
        Imports data from apex into own db
        """
        self.statusBar().showMessage('')
        self.close_db()
        self.start_db()
        self.ui.cifList_treeWidget.show()
        self.abort_import_button.show()
        n = 1
        num = 0
        time1 = time.perf_counter()
        conn = self.open_apex_db(user, password, host)
        if not conn:
            self.abort_import_button.hide()
            return None
        cif = Cif()
        with self.session_scope() as session:
            for apd in self.apx.get_all_data():
                if num == 20:
                    num = 0
                self.progressbar(num, 0, 20)
                cif.cif_data['_cell_length_a'] = apd[1]
                cif.cif_data['_cell_length_b'] = apd[2]
                cif.cif_data['_cell_length_c'] = apd[3]
                cif.cif_data['_cell_angle_alpha'] = apd[4]
                cif.cif_data['_cell_angle_beta'] = apd[5]
                cif.cif_data['_cell_angle_gamma'] = apd[6]
                cif.cif_data["data"] = apd[8]
                cif.cif_data['_diffrn_radiation_wavelength'] = apd[13]
                cif.cif_data['_exptl_crystal_colour'] = apd[14]
                cif.cif_data['_exptl_crystal_size_max'] = apd[16]
                cif.cif_data['_exptl_crystal_size_mid'] = apd[17]
                cif.cif_data['_exptl_crystal_size_min'] = apd[18]
                cif.cif_data["_chemical_formula_sum"] = apd[25]
                cif.cif_data['_diffrn_reflns_av_R_equivalents'] = apd[21]  # rint
                cif.cif_data['_diffrn_reflns_av_unetI/netI'] = apd[22]  # rsig
                cif.cif_data['_diffrn_reflns_number'] = apd[23]
                comp = apd[26]
                if comp:
                    cif.cif_data['_diffrn_measured_fraction_theta_max'] = comp / 100
                tst = fill_db_with_cif_data(session=session, cif=cif, filename=apd[8], path=apd[12], structure_id=n)
                if not tst:
                    continue
                self.add_table_row(filename=apd[8], data=apd[8], path=apd[12], structure_id=str(n))
                n += 1
                if n % 300 == 0:
                    pass
                    # session.commit()
                num += 1
                if not self.decide_import:
                    # This means, import was aborted.
                    self.abort_import_button.hide()
                    self.decide_import = True
                    break
        time2 = time.perf_counter()
        diff = time2 - time1
        self.progress.hide()
        m, s = divmod(diff, 60)
        h, m = divmod(m, 60)
        self.ui.statusbar.showMessage('Added {} APEX entries in: {:>2d} h, {:>2d} m, {:>3.2f} s'
                                      .format(n, int(h), int(m), s), msecs=0)
        # self.ui.cifList_treeWidget.resizeColumnToContents(0)
        self.set_columnsize()
        init_textsearch(self.engine)
        populate_fulltext_search_table(self.engine)
        self.abort_import_button.hide()

    def set_columnsize(self):
        """
        Sets columnsize of main structure list.
        """
        self.ui.cifList_treeWidget.sortByColumn(0, 0)
        treewidth = self.ui.cifList_treeWidget.width()
        self.ui.cifList_treeWidget.setColumnWidth(0, treewidth / 4)
        self.ui.cifList_treeWidget.setColumnWidth(1, treewidth / 5)
        # self.ui.cifList_treeWidget.resizeColumnToContents(0)
        # self.ui.cifList_treeWidget.resizeColumnToContents(1)

    def show_full_list(self) -> None:
        """
        Displays the complete list of structures
        [structure_id, meas, path, filename, data]
        """
        self.ui.cifList_treeWidget.clear()
        las_rowid = 0
        with self.session_scope() as session:
            las_rowid = get_lastrow_id(session)
            for Id, filename, dataname, path in session.query(Structure.Id, Structure.filename,
                                                              Structure.dataname, Structure.path).all():
                self.add_table_row(filename, path, dataname, Id)
        mess = "Loaded {} entries.".format(las_rowid)
        self.statusBar().showMessage(mess, msecs=5000)
        self.set_columnsize()
        # self.ui.cifList_treeWidget.resizeColumnToContents(0)
        # self.ui.cifList_treeWidget.resizeColumnToContents(1)
        self.full_list = True
        self.ui.tabWidget.setCurrentIndex(0)
        self.ui.SpGrcomboBox.setCurrentIndex(0)
        self.ui.ad_elementsIncLineEdit.clear()
        self.ui.ad_elementsExclLineEdit.clear()
        self.ui.ad_moreResultscheckBox.setChecked(False)
        self.ui.ad_superlatticeCheckBox.setChecked(False)
        self.ui.ad_textsearch.clear()
        self.ui.ad_textsearch_excl.clear()
        self.ui.ad_unitCellLineEdit.clear()
        self.ui.dateEdit1.setDate(QtCore.QDate(date.today()))
        self.ui.dateEdit2.setDate(QtCore.QDate(date.today()))

    def clear_fields(self) -> None:
        """
        Clears all residuals fields.
        """
        self.ui.completeLineEdit.clear()
        self.ui.dataReflnsLineEdit.clear()
        self.ui.dLineEdit.clear()
        self.ui.goofLineEdit.clear()
        self.ui.maxShiftLineEdit.clear()
        self.ui.numParametersLineEdit.clear()
        self.ui.numRestraintsLineEdit.clear()
        self.ui.peakLineEdit.clear()
        self.ui.r1LineEdit.clear()
        self.ui.reflTotalLineEdit.clear()
        self.ui.rintLineEdit.clear()
        self.ui.rsigmaLineEdit.clear()
        self.ui.SpaceGroupLineEdit.clear()
        self.ui.formLabel.clear()
        self.ui.temperatureLineEdit.clear()
        self.ui.thetaFullLineEdit.clear()
        self.ui.thetaMaxLineEdit.clear()
        self.ui.wavelengthLineEdit.clear()
        self.ui.wR2LineEdit.clear()
        self.ui.zLineEdit.clear()
        self.ui.cCDCNumberLineEdit.clear()
        self.ui.refl2sigmaLineEdit.clear()
        self.ui.uniqReflLineEdit.clear()
        self.ui.lastModifiedLineEdit.clear()


if __name__ == "__main__":
    try:
        uic.compileUiDir(os.path.join(application_path, './gui'))
    except:
        print("Unable to compile UI!")
    from gui.strf_main import Ui_stdbMainwindow
    from gui.strf_dbpasswd import Ui_PasswdDialog

    # later http://www.pyinstaller.org/
    app = QtWidgets.QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon('./icons/strf.png'))
    # Has to be without version number, because QWebengine stores data in ApplicationName directory:
    app.setApplicationName('StructureFinder')
    # app.setApplicationDisplayName("StructureFinder")
    myapp = StartStructureDB()
    myapp.show()
    myapp.raise_()
    myapp.setWindowTitle('StructureFinder v{}'.format(VERSION))
    try:
        sys.exit(app.exec_())
    except KeyboardInterrupt:
        sys.exit()
