from __future__ import print_function

import os
from pprint import pprint
import sys

import time
from PyQt5 import uic, Qt3DExtras, QtWidgets, Qt3DRender
from PyQt5.QtCore import QSize, pyqtSignal, pyqtSlot
from PyQt5.QtGui import QColor, QVector3D
from PyQt5.QtQml import QQmlApplicationEngine
from PyQt5.QtWidgets import QApplication, QWidget
from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QTreeWidgetItem

from lattice import lattice
from searcher import filecrawler
from stdb_main import Ui_stdbMainwindow
from searcher.fileparser import Cif
from searcher.database_handler import StructureTable, DatabaseRequest
from opengl.moleculegl import MyScene

uic.compileUiDir('./')


"""
TODO:
- make progress bar for indexer and file opener
- structure code
- make 3D model from atoms
- make file type more flexible. handle .res and .cif equally
- group structures in measurements
- list properties of a selected cif file
- implement progress bar for indexing
- implement "save on close?" dialog
- add abort button for indexer
- recognize already indexed files
- search for strings to get a result for a persons name, add person to db
"""


class StartStructureDB(QMainWindow):
    #changedValue = pyqtSignal('QString')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ui = Ui_stdbMainwindow()
        self.ui.setupUi(self)
        self.connect_signals_and_slots()
        self.statusBar().showMessage('Ready')
        self.ui.cifList_treeWidget.show()
        self.ui.cifList_treeWidget.hideColumn(2)
        # self.ui.cellSearchEdit.hide()
        self.dbfilename = 'test.sqlite'
        print(self.dbfilename)

        #######################################################
        """
        # TODO: pull this out:
        view = Qt3DExtras.Qt3DWindow()
        s = MyScene()
        scene = s.createScene()
        print('#scene')
        # // Camera
        camera = view.camera()
        lens = Qt3DRender.QCameraLens()
        lens.setPerspectiveProjection(45.0, 16.0 / 9.0, 0.1, 1000.0)
        camera.setProjectionType(Qt3DRender.QCameraLens.PerspectiveProjection)
        camera.setUpVector(QVector3D(0, 1.0, 0))
        camera.setPosition(QVector3D(0, 0, 140.0))  # Entfernung
        camera.setViewCenter(QVector3D(0, 0, 0))
        print('#camera')
        # // For camera controls
        camController = Qt3DExtras.QOrbitCameraController(scene)
        camController.setLinearSpeed(50.0)
        camController.setLookSpeed(180.0)
        camController.setCamera(camera)
        #view.setRootEntity(scene)
        print('view#')
        view.defaultFrameGraph().setClearColor(QColor('lightgray'))
        container = QWidget.createWindowContainer(view)
        screenSize = view.screen().size()
        container.setMinimumSize(QSize(200, 100))
        container.setMaximumSize(screenSize)
        self.ui.openglVlayout.addWidget(container, 1)
        view.show()
        #########################################
        """
        self.ui.centralwidget.setMinimumSize(1200, 500)
        self.showMaximized()
        try:
            # TODO: don't do in future:
            os.remove(self.dbfilename)
        except:
            pass
        self.structures = StructureTable(self.dbfilename)
        self.db = DatabaseRequest(self.dbfilename)
        self.db.initialize_db()
        self.show()

    def connect_signals_and_slots(self):
        self.ui.importDatabaseButton.clicked.connect(self.import_database)
        self.ui.importDirButton.clicked.connect(self.import_cif_dirs)
        self.ui.searchLineEDit.textChanged.connect(self.search_cell)
        # self.ui.actionExit.triggered.connect(QtGui.QGuiApplication.quit)
        self.ui.cifList_treeWidget.clicked.connect(self.show_properties)
        # for later use to implement relocation of whole database:
        # das brauch ich nicht:
        # self.ui.cifList_treeWidget.doubleClicked.connect(self.relocate)
        # self.ui.cifList_treeWidget.doubleClicked.connect(self.show_properties)

    def show_properties(self, item):
        """
        This slot show the properties of a cif file in the properties widget
        """
        # self.ui.properties_treeWidget.show()
        cell = self.structures.get_cell_by_id(item.sibling(item.row(), 2).data())
        # print(item.sibling(item.row(), 2).data())
        a, b, c, alpha, beta, gamma = cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]
        if a:
            self.ui.aLineEdit.setText("{:>5.4f}".format(a))
        if b:
            self.ui.bLineEdit.setText("{:>5.4f}".format(b))
        if c:
            self.ui.cLineEdit.setText("{:>5.4f}".format(c))
        if alpha:
            self.ui.alphaLineEdit.setText("{:>5.4f}".format(alpha))
        if beta:
            self.ui.betaLineEdit.setText("{:>5.4f}".format(beta))
        if gamma:
            self.ui.gammaLineEdit.setText("{:>5.4f}".format(gamma))

    @pyqtSlot('QString')
    def search_cell(self, search_string):
        """
        searches db for given cell via the cell volume
        
        8.4009  10.4848  11.8979  94.7910 103.0250 108.5480
        
        :param search_string: 
        :return: 
        """
        try:
            cell = [float(x) for x in search_string.split()]
        except (TypeError, ValueError):
            return False
        if len(cell) != 6:
            return False
        try:
            self.structures = StructureTable(self.dbfilename)
            volume = lattice.vol_unitcell(*cell)
            res = self.structures.find_by_volume(volume)
            print(res)
        except ValueError:
            return False

    def import_database(self):
        """
        Import a new database.
        :return: 
        """
        fname = QFileDialog.getOpenFileName(self, 'Open File', '')
        print(fname)
        self.dbfilename = fname[0]
        self.structures = StructureTable(self.dbfilename)
        self.ui.cifList_treeWidget.show()
        if not self.structures:
            return False
        for i in self.structures.get_all_structure_names():
            """structure.Id, structure.measurement, structure.path, structure.filename, 
                         structure.dataname"""
            str_tree = QTreeWidgetItem(self.ui.cifList_treeWidget)
            str_tree.setText(0, i[3])  # name
            str_tree.setText(1, i[2])  # path
            str_tree.setData(2, 0, i[0])  # id
        self.ui.cifList_treeWidget.resizeColumnToContents(0)

    def import_cif_dirs(self):
        """
        Imports cif files from a certain directory
        :return: None
        """
        fname = QFileDialog.getExistingDirectory(self, 'Open Directory', '')
        #fname = "/Users/daniel/Documents/Strukturen/Miriam/IKms_cf_08_Ni(mes)(cod)PF/FINAL/"
        # fname = "D:/GitHub/StructureDB/test-data"
        # fname = os.path.abspath("/Users/daniel/Downloads")
        # fname = os.path.abspath("../")
        if not fname:
            return False
        time1 = time.clock()
        try:
            files = list(filecrawler.create_file_list(str(fname), endings='cif'))
        except FileNotFoundError as e:
            print(e)
            return False
        time2 = time.clock()
        diff = time2 - time1
        print("File list:", round(diff, 4), 's')
        self.ui.cifList_treeWidget.show()
        # TODO: implement multiple cells in one cif file:
        n = 1
        times = []
        for filepth in files:
            filename = filepth.name
            path = str(filepth.parents[0])
            #print(filepth)  # print full file path
            structure_id = n
            time2 = time.clock()
            cif = Cif(filepth)
            time3 = time.clock()
            diff2 = time3 - time2
            times.append(diff2)
            #print(round(diff, 4), 's')
            if not cif.ok:
                continue
            #print(cif, '##')
            if cif and filename and path:
                self.fill_tables(cif, filename, path, structure_id)
                strTree = QTreeWidgetItem(self.ui.cifList_treeWidget)
                strTree.setText(0, filename)
                strTree.setText(1, path)
                strTree.setText(2, str(n))
                n += 1
        print('Parse cif files:', round(sum(times), 3), 's,', n, 'files')
        self.ui.cifList_treeWidget.resizeColumnToContents(0)
        #self.ui.cifList_treeWidget.resizeColumnToContents(1)
        # self.ui.relocate_lineEdit.hide()
        self.structures.database.commit_db("Committed")

    def fill_tables(self, cif, filename, path, structure_id):
        """
        FIll all info from cif file into the database tables 
        :param cif: 
        :param filename: 
        :param path: 
        :param structure_id: 
        :return: 
        """
        a = cif._cell_length_a
        b = cif._cell_length_b
        c = cif._cell_length_c
        alpha = cif._cell_angle_alpha
        beta = cif._cell_angle_beta
        gamma = cif._cell_angle_gamma
        measurement_id = self.structures.fill_measuremnts_table(filename, structure_id)
        self.structures.fill_structures_table(path, filename, structure_id, measurement_id, cif.cif_data['data'])
        self.structures.fill_cell_table(structure_id, a, b, c, alpha, beta, gamma)
        #pprint(cif._atom)
        for x in cif._atom:
            try:
                self.structures.fill_atoms_table(structure_id, x,
                                             cif._atom[x]['_atom_site_type_symbol'],
                                             cif._atom[x]['_atom_site_fract_x'].split('(')[0],
                                             cif._atom[x]['_atom_site_fract_y'].split('(')[0],
                                             cif._atom[x]['_atom_site_fract_z'].split('(')[0])
            except KeyError as e:
                pass
                #print("Atom:", x, path, filename)
                #print(e)



class QmlAusgabe(object):
    def __init__(self, pathToQmlFile="beispiel.qml"):
        #QML-Engine
        self.__appEngine = QQmlApplicationEngine()
        self.__appEngine.load(pathToQmlFile)
        self.__appWindow = self.__appEngine.rootObjects()[0]

    def show(self):
        self.__appWindow.show()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    myapp = StartStructureDB()
    myapp.show()
    myapp.raise_()
    try:
        sys.exit(app.exec_())
    except KeyboardInterrupt:
        sys.exit()
