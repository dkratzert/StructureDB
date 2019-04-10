"""
Unit tests for StructureFinder
"""
import doctest
import platform
import sys
import unittest
from contextlib import suppress
from pathlib import Path

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon
from PyQt5.QtTest import QTest
from PyQt5.QtWidgets import QApplication

import searcher
import searcher.misc
import strf
from misc import update_check
from misc.version import VERSION
from pymatgen.core import lattice
from searcher import database_handler, fileparser
from shelxfile import shelx, elements, misc


class DoctestsTest(unittest.TestCase):
    def testrun_doctest(self):
        for name in [strf, shelx, elements, misc, searcher, update_check, database_handler,
                     fileparser, searcher.misc]:
            failed, attempted = doctest.testmod(name)  # , verbose=True)
            if failed == 0:
                print('passed all {} tests in {}!'.format(attempted, name.__name__))
            else:
                msg = '!!!!!!!!!!!!!!!! {} of {} tests failed in {}  !!!!!!!!!!!!!!!!!!!!!!!!!!!'.format(failed,
                                                                                                         attempted,
                                                                                                         name.__name__)
                self.assertFalse(failed, msg)


app = QApplication(sys.argv)


class TestApplication(unittest.TestCase):

    def setUp(self) -> None:
        # uic.compileUiDir('./gui')
        app.setWindowIcon(QIcon('./icons/strf.png'))
        # Has to be without version number, because QWebengine stores data in ApplicationName directory:
        app.setApplicationName('StructureFinder')
        self.myapp = strf.StartStructureDB()
        self.myapp.setWindowTitle('StructureFinder v{}'.format(VERSION))
        self.myapp.structures = database_handler.StructureTable('./test-data/test.sql')
        self.myapp.show_full_list()

    def tearDown(self) -> None:
        super(TestApplication, self).tearDown()

    # @unittest.skip("foo")
    def test_gui_simpl(self):
        # Number of items in main list
        self.assertEqual(263, self.myapp.ui.cifList_treeWidget.topLevelItemCount())
        # structureId
        self.assertEqual('241', self.myapp.ui.cifList_treeWidget.topLevelItem(1).text(3))
        # filename
        self.assertEqual('1000000.cif', self.myapp.ui.cifList_treeWidget.topLevelItem(1).text(0))

    # @unittest.skip('skipping unfinished')
    def test_search_cell_simpl(self):
        """
        Testing simple unit cell search.
        """
        # correct cell:
        self.myapp.ui.searchCellLineEDit.setText('7.878 10.469 16.068 90.000 95.147 90.000')
        self.assertEqual(3, self.myapp.ui.cifList_treeWidget.topLevelItemCount())
        self.myapp.show_full_list()
        # incomplete unit cell:
        self.myapp.ui.searchCellLineEDit.setText('7.878 10.469 16.068 90.000 95.147')
        self.assertEqual(263, self.myapp.ui.cifList_treeWidget.topLevelItemCount())
        self.myapp.show_full_list()
        # invalid unit cell:
        self.myapp.ui.searchCellLineEDit.setText('7.878 10.469 16.068 90.000 95.147 abc')
        self.assertEqual(263, self.myapp.ui.cifList_treeWidget.topLevelItemCount())
        self.assertEqual("Not a valid unit cell!", self.myapp.statusBar().currentMessage())

    def test_search_text_simpl(self):
        """
        Testing simple text search.
        """
        self.myapp.ui.txtSearchEdit.setText('SADI')
        self.assertEqual(4, self.myapp.ui.cifList_treeWidget.topLevelItemCount())
        self.assertEqual("Found 4 entries.", self.myapp.statusBar().currentMessage())
        self.myapp.show_full_list()
        self.myapp.ui.txtSearchEdit.setText('sadi')
        self.assertEqual(4, self.myapp.ui.cifList_treeWidget.topLevelItemCount())
        self.myapp.show_full_list()
        # should give no result
        self.myapp.ui.txtSearchEdit.setText('foobar')
        self.assertEqual(0, self.myapp.ui.cifList_treeWidget.topLevelItemCount())
        self.assertEqual("Found 0 entries.", self.myapp.statusBar().currentMessage())

    # @unittest.skip("foo")
    def test_clicks(self):
        """
        Testing copy to clip board with double click on unit cell
        """
        item = self.myapp.ui.cifList_treeWidget.topLevelItem(0)
        self.myapp.ui.cifList_treeWidget.setCurrentItem(item)
        QTest.mouseDClick(self.myapp.ui.cellField, Qt.LeftButton, delay=5)
        clp = QApplication.clipboard().text()
        self.assertEqual(" 7.878 10.469 16.068 90.000 95.147 90.000", clp)

    def test_save_db(self):
        """
        Saves the current database to a file.
        """
        self.myapp.import_file_dirs('test-data/COD')
        testfile = Path('./tst.sql')
        with suppress(Exception):
            Path.unlink(testfile)
        self.myapp.save_database(testfile.absolute())
        self.assertEqual(True, testfile.is_file())
        self.assertEqual(True, testfile.exists())
        Path.unlink(testfile)
        self.assertEqual(False, testfile.exists())
        self.assertEqual('Database saved.', self.myapp.statusBar().currentMessage())

    def test_index_db(self):
        """
        Test index and save
        """
        self.myapp.import_file_dirs('test-data/COD')
        self.assertEqual(22, self.myapp.ui.cifList_treeWidget.topLevelItemCount())
        self.myapp.import_file_dirs('gui')
        self.assertEqual(0, self.myapp.ui.cifList_treeWidget.topLevelItemCount())
        self.myapp.import_file_dirs('test-data/tst')
        self.assertEqual(3, self.myapp.ui.cifList_treeWidget.topLevelItemCount())
        self.myapp.ui.add_cif.setChecked(False)
        self.myapp.ui.add_res.setChecked(True)
        self.myapp.import_file_dirs('test-data/tst')
        self.assertEqual(1, self.myapp.ui.cifList_treeWidget.topLevelItemCount())

    def test_open_database_file(self):
        """
        Testing the opening of a database.
        """
        # self.myapp.close_db()  # not needed here!
        status = self.myapp.import_database_file('test-data/test.sql')
        self.assertEqual(True, status)
        self.assertEqual(263, self.myapp.ui.cifList_treeWidget.topLevelItemCount())

    def test_p4p_parser(self):
        self.myapp.search_for_p4pcell('test-data/test2.p4p')
        self.assertEqual('14.637 9.221  15.094 90.000 107.186 90.000', self.myapp.ui.searchCellLineEDit.text())
        # self.assertEqual(263, self.myapp.ui.cifList_treeWidget.topLevelItemCount())

    def test_res_parser(self):
        self.myapp.search_for_res_cell('test-data/p21c.res')
        self.assertEqual('10.509 20.904 20.507 90.000 94.130 90.000', self.myapp.ui.searchCellLineEDit.text())
        self.assertEqual(2, self.myapp.ui.cifList_treeWidget.topLevelItemCount())

    def test_all_cif_values(self):
        item = self.myapp.ui.cifList_treeWidget.topLevelItem(3)
        self.myapp.ui.cifList_treeWidget.setCurrentItem(item)
        QTest.mouseClick(self.myapp.ui.allEntrysTab, Qt.LeftButton)
        self.assertEqual('Id', self.myapp.ui.allCifTreeWidget.topLevelItem(0).text(0))
        self.assertEqual('250', self.myapp.ui.allCifTreeWidget.topLevelItem(0).text(1))
        self.assertEqual('C107 H142 N14 O26', self.myapp.ui.allCifTreeWidget.topLevelItem(10).text(1))

    def test_res_file_tab(self):
        item = self.myapp.ui.cifList_treeWidget.topLevelItem(262)
        self.myapp.ui.cifList_treeWidget.setCurrentItem(item)
        QTest.mouseClick(self.myapp.ui.SHELXtab, Qt.LeftButton)
        self.assertEqual('REM Solution', self.myapp.ui.SHELXplainTextEdit.toPlainText()[:12])
        ###############
        item = self.myapp.ui.cifList_treeWidget.topLevelItem(260)
        self.myapp.ui.cifList_treeWidget.setCurrentItem(item)
        QTest.mouseClick(self.myapp.ui.SHELXtab, Qt.LeftButton)
        self.assertEqual('TITL p21c in', self.myapp.ui.SHELXplainTextEdit.toPlainText()[:12])
        ###############
        item = self.myapp.ui.cifList_treeWidget.topLevelItem(250)
        self.myapp.ui.cifList_treeWidget.setCurrentItem(item)
        QTest.mouseClick(self.myapp.ui.SHELXtab, Qt.LeftButton)
        self.assertEqual('No SHELXL res file in cif found.', self.myapp.ui.SHELXplainTextEdit.toPlainText())

    def test_cellchekcsd(self):
        item = self.myapp.ui.cifList_treeWidget.topLevelItem(2)
        self.myapp.ui.cifList_treeWidget.setCurrentItem(item)
        QTest.mouseClick(self.myapp.ui.CCDCSearchTab, Qt.LeftButton)
        if platform.system() == 'Windows':
            self.assertEqual(' 7.878 10.469 16.068 90.000 95.147 90.000', self.myapp.ui.cellSearchLabelCSD.text())

    def test_adv(self):
        pass

######################################################
##  Database testing:
######################################################


class TestSearch(unittest.TestCase):
    def setUp(self) -> None:
        self.dbfilename = 'test-data/test.sql'
        self.structures = database_handler.StructureTable(self.dbfilename)
        # more results:
        self.m_vol_threshold = 0.04
        self.m_ltol = 0.08
        self.m_atol = 1.8
        # regular:
        self.vol_threshold = 0.02
        self.ltol = 0.03
        self.atol = 1.0

    def test_cellfind(self):
        idlist = []
        cell = [10.930, 12.716, 15.709, 90.000, 90.000, 90.000]
        results = [(49, 9.242, 12.304, 18.954, 90.0, 92.74, 90.0, 2153.0),
                   (260, 10.93, 12.7162, 15.7085, 90.0, 90.0, 90.0, 2183.3),
                   (10, 13.5918, 10.7345, 16.442, 90.0, 113.142, 90.0, 2205.9),
                   (244, 13.5918, 10.7345, 16.442, 90.0, 113.142, 90.0, 2205.9),
                   (207, 16.139, 5.117, 26.887, 90.0, 90.0, 90.0, 2220.4)]
        volume = searcher.misc.vol_unitcell(*cell)
        cells = self.structures.find_by_volume(volume, self.vol_threshold)
        self.assertEqual(cells, results)
        lattice1 = lattice.Lattice.from_parameters_niggli_reduced(*cell)
        for curr_cell in cells:
            try:
                lattice2 = lattice.Lattice.from_parameters(*curr_cell[1:7])
            except ValueError:
                continue
            mapping = lattice1.find_mapping(lattice2, self.ltol, self.atol, skip_rotation_matrix=True)
            if mapping:
                idlist.append(curr_cell[0])
        self.assertEqual(idlist, [260])

    def test_more_results_cellfind(self):
        idlist = []
        cell = [10.930, 12.716, 15.709, 90.000, 90.000, 90.000]
        results = [(251, 13.432, 10.5988, 16.2393, 90.0, 113.411, 90.0, 2121.6),
                   (161, 14.8208, 8.1939, 17.4844, 90.0, 91.185, 90.0, 2122.9),
                   (49, 9.242, 12.304, 18.954, 90.0, 92.74, 90.0, 2153.0),
                   (260, 10.93, 12.7162, 15.7085, 90.0, 90.0, 90.0, 2183.3),
                   (10, 13.5918, 10.7345, 16.442, 90.0, 113.142, 90.0, 2205.9),
                   (244, 13.5918, 10.7345, 16.442, 90.0, 113.142, 90.0, 2205.9),
                   (207, 16.139, 5.117, 26.887, 90.0, 90.0, 90.0, 2220.4),
                   (71, 14.815, 14.264, 10.55, 90.0, 90.0, 90.0, 2229.4),
                   (113, 15.187, 12.883, 11.468, 90.0, 90.0, 90.0, 2243.8),
                   (129, 27.858, 8.094, 9.951, 90.0, 90.0, 90.0, 2243.8),
                   (1, 10.36, 18.037, 25.764, 127.03, 129.81, 90.51, 2260.487670154818),
                   (12, 10.36, 18.037, 25.764, 127.03, 129.81, 90.51, 2260.487670154818)]
        volume = searcher.misc.vol_unitcell(*cell)
        cells = self.structures.find_by_volume(volume, self.m_vol_threshold)
        self.assertEqual(cells, results)
        lattice1 = lattice.Lattice.from_parameters_niggli_reduced(*cell)
        for curr_cell in cells:
            try:
                lattice2 = lattice.Lattice.from_parameters(*curr_cell[1:7])
            except ValueError:
                continue
            mapping = lattice1.find_mapping(lattice2, self.m_ltol, self.m_atol, skip_rotation_matrix=True)
            if mapping:
                idlist.append(curr_cell[0])
        self.assertEqual(idlist, [260, 113])
