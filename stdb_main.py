# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file './stdb_main.ui'
#
# Created by: PyQt5 UI code generator 5.8.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_stdbMainwindow(object):
    def setupUi(self, stdbMainwindow):
        stdbMainwindow.setObjectName("stdbMainwindow")
        stdbMainwindow.resize(958, 659)
        self.centralwidget = QtWidgets.QWidget(stdbMainwindow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.dataDock = QtWidgets.QDockWidget(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.dataDock.sizePolicy().hasHeightForWidth())
        self.dataDock.setSizePolicy(sizePolicy)
        self.dataDock.setMinimumSize(QtCore.QSize(400, 345))
        self.dataDock.setFloating(False)
        self.dataDock.setFeatures(QtWidgets.QDockWidget.DockWidgetMovable)
        self.dataDock.setObjectName("dataDock")
        self.dockWidgetContents_3 = QtWidgets.QWidget()
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(6)
        sizePolicy.setVerticalStretch(6)
        sizePolicy.setHeightForWidth(self.dockWidgetContents_3.sizePolicy().hasHeightForWidth())
        self.dockWidgetContents_3.setSizePolicy(sizePolicy)
        self.dockWidgetContents_3.setMinimumSize(QtCore.QSize(100, 100))
        self.dockWidgetContents_3.setObjectName("dockWidgetContents_3")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.dockWidgetContents_3)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 13)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.groupBox = QtWidgets.QGroupBox(self.dockWidgetContents_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        self.groupBox.setTitle("")
        self.groupBox.setObjectName("groupBox")
        self.gridLayout = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout.setObjectName("gridLayout")
        self.importDatabaseButton = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.importDatabaseButton.sizePolicy().hasHeightForWidth())
        self.importDatabaseButton.setSizePolicy(sizePolicy)
        self.importDatabaseButton.setObjectName("importDatabaseButton")
        self.gridLayout.addWidget(self.importDatabaseButton, 2, 0, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem, 3, 0, 1, 1)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem1, 10, 0, 1, 1)
        spacerItem2 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem2, 7, 0, 1, 1)
        self.optionsButton = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.optionsButton.sizePolicy().hasHeightForWidth())
        self.optionsButton.setSizePolicy(sizePolicy)
        self.optionsButton.setObjectName("optionsButton")
        self.gridLayout.addWidget(self.optionsButton, 11, 0, 1, 1)
        self.searchButton = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.searchButton.sizePolicy().hasHeightForWidth())
        self.searchButton.setSizePolicy(sizePolicy)
        self.searchButton.setObjectName("searchButton")
        self.gridLayout.addWidget(self.searchButton, 4, 0, 1, 1)
        self.importDirButton = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.importDirButton.sizePolicy().hasHeightForWidth())
        self.importDirButton.setSizePolicy(sizePolicy)
        self.importDirButton.setObjectName("importDirButton")
        self.gridLayout.addWidget(self.importDirButton, 0, 0, 1, 1)
        spacerItem3 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem3, 1, 0, 1, 1)
        self.exportButton = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.exportButton.sizePolicy().hasHeightForWidth())
        self.exportButton.setSizePolicy(sizePolicy)
        self.exportButton.setObjectName("exportButton")
        self.gridLayout.addWidget(self.exportButton, 8, 0, 1, 1)
        self.gridLayout_2.addWidget(self.groupBox, 1, 0, 1, 1)
        self.stackedWidget = QtWidgets.QStackedWidget(self.dockWidgetContents_3)
        self.stackedWidget.setObjectName("stackedWidget")
        self.importResults = QtWidgets.QWidget()
        self.importResults.setObjectName("importResults")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.importResults)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.relocate_lineEdit = QtWidgets.QLineEdit(self.importResults)
        self.relocate_lineEdit.setFrame(False)
        self.relocate_lineEdit.setObjectName("relocate_lineEdit")
        self.verticalLayout.addWidget(self.relocate_lineEdit)
        self.cifList_treeWidget = QtWidgets.QTreeWidget(self.importResults)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(6)
        sizePolicy.setVerticalStretch(4)
        sizePolicy.setHeightForWidth(self.cifList_treeWidget.sizePolicy().hasHeightForWidth())
        self.cifList_treeWidget.setSizePolicy(sizePolicy)
        self.cifList_treeWidget.setMinimumSize(QtCore.QSize(0, 288))
        self.cifList_treeWidget.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.cifList_treeWidget.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.cifList_treeWidget.setObjectName("cifList_treeWidget")
        self.verticalLayout.addWidget(self.cifList_treeWidget)
        self.properties_treeWidget = QtWidgets.QTreeWidget(self.importResults)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(5)
        sizePolicy.setHeightForWidth(self.properties_treeWidget.sizePolicy().hasHeightForWidth())
        self.properties_treeWidget.setSizePolicy(sizePolicy)
        self.properties_treeWidget.setObjectName("properties_treeWidget")
        item_0 = QtWidgets.QTreeWidgetItem(self.properties_treeWidget)
        item_0 = QtWidgets.QTreeWidgetItem(self.properties_treeWidget)
        item_1 = QtWidgets.QTreeWidgetItem(item_0)
        item_1 = QtWidgets.QTreeWidgetItem(item_0)
        item_1 = QtWidgets.QTreeWidgetItem(item_0)
        item_1 = QtWidgets.QTreeWidgetItem(item_0)
        item_1 = QtWidgets.QTreeWidgetItem(item_0)
        item_1 = QtWidgets.QTreeWidgetItem(item_0)
        item_0 = QtWidgets.QTreeWidgetItem(self.properties_treeWidget)
        item_0 = QtWidgets.QTreeWidgetItem(self.properties_treeWidget)
        item_0 = QtWidgets.QTreeWidgetItem(self.properties_treeWidget)
        item_0 = QtWidgets.QTreeWidgetItem(self.properties_treeWidget)
        item_1 = QtWidgets.QTreeWidgetItem(item_0)
        item_1 = QtWidgets.QTreeWidgetItem(item_0)
        item_0 = QtWidgets.QTreeWidgetItem(self.properties_treeWidget)
        item_1 = QtWidgets.QTreeWidgetItem(item_0)
        item_1 = QtWidgets.QTreeWidgetItem(item_0)
        item_0 = QtWidgets.QTreeWidgetItem(self.properties_treeWidget)
        item_0 = QtWidgets.QTreeWidgetItem(self.properties_treeWidget)
        item_0 = QtWidgets.QTreeWidgetItem(self.properties_treeWidget)
        item_0 = QtWidgets.QTreeWidgetItem(self.properties_treeWidget)
        self.verticalLayout.addWidget(self.properties_treeWidget)
        self.stackedWidget.addWidget(self.importResults)
        self.refinement_edit_page = QtWidgets.QWidget()
        self.refinement_edit_page.setObjectName("refinement_edit_page")
        self.moleculeGLWidget = QtWidgets.QOpenGLWidget(self.refinement_edit_page)
        self.moleculeGLWidget.setEnabled(False)
        self.moleculeGLWidget.setGeometry(QtCore.QRect(0, 0, 721, 541))
        self.moleculeGLWidget.setObjectName("moleculeGLWidget")
        self.stackedWidget.addWidget(self.refinement_edit_page)
        self.stackedWidgetPage2 = QtWidgets.QWidget()
        self.stackedWidgetPage2.setEnabled(False)
        self.stackedWidgetPage2.setObjectName("stackedWidgetPage2")
        self.stackedWidget.addWidget(self.stackedWidgetPage2)
        self.gridLayout_2.addWidget(self.stackedWidget, 1, 1, 1, 1)
        self.dataDock.setWidget(self.dockWidgetContents_3)
        self.horizontalLayout.addWidget(self.dataDock)
        stdbMainwindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(stdbMainwindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 958, 22))
        self.menubar.setObjectName("menubar")
        self.fileMenu = QtWidgets.QMenu(self.menubar)
        self.fileMenu.setObjectName("fileMenu")
        self.menuEdit = QtWidgets.QMenu(self.menubar)
        self.menuEdit.setObjectName("menuEdit")
        self.menuHelp = QtWidgets.QMenu(self.menubar)
        self.menuHelp.setObjectName("menuHelp")
        stdbMainwindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(stdbMainwindow)
        self.statusbar.setObjectName("statusbar")
        stdbMainwindow.setStatusBar(self.statusbar)
        self.actionExit = QtWidgets.QAction(stdbMainwindow)
        self.actionExit.setObjectName("actionExit")
        self.actionEdit_dataset = QtWidgets.QAction(stdbMainwindow)
        self.actionEdit_dataset.setObjectName("actionEdit_dataset")
        self.actionUser_manual = QtWidgets.QAction(stdbMainwindow)
        self.actionUser_manual.setObjectName("actionUser_manual")
        self.actionImport_file = QtWidgets.QAction(stdbMainwindow)
        self.actionImport_file.setObjectName("actionImport_file")
        self.actionImport_directory = QtWidgets.QAction(stdbMainwindow)
        self.actionImport_directory.setObjectName("actionImport_directory")
        self.actionExport_Database_s = QtWidgets.QAction(stdbMainwindow)
        self.actionExport_Database_s.setObjectName("actionExport_Database_s")
        self.actionOptions = QtWidgets.QAction(stdbMainwindow)
        self.actionOptions.setObjectName("actionOptions")
        self.fileMenu.addAction(self.actionExit)
        self.fileMenu.addAction(self.actionImport_file)
        self.fileMenu.addAction(self.actionImport_directory)
        self.fileMenu.addAction(self.actionExport_Database_s)
        self.fileMenu.addAction(self.actionOptions)
        self.menuEdit.addAction(self.actionEdit_dataset)
        self.menuHelp.addAction(self.actionUser_manual)
        self.menubar.addAction(self.fileMenu.menuAction())
        self.menubar.addAction(self.menuEdit.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(stdbMainwindow)
        self.stackedWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(stdbMainwindow)

    def retranslateUi(self, stdbMainwindow):
        _translate = QtCore.QCoreApplication.translate
        stdbMainwindow.setWindowTitle(_translate("stdbMainwindow", "MainWindow"))
        self.importDatabaseButton.setText(_translate("stdbMainwindow", "Open Database"))
        self.optionsButton.setText(_translate("stdbMainwindow", "Options"))
        self.searchButton.setText(_translate("stdbMainwindow", "Cell Search"))
        self.importDirButton.setText(_translate("stdbMainwindow", "Import Directory"))
        self.exportButton.setText(_translate("stdbMainwindow", "Export Database(s)"))
        self.cifList_treeWidget.setSortingEnabled(True)
        self.cifList_treeWidget.headerItem().setText(0, _translate("stdbMainwindow", "file"))
        self.cifList_treeWidget.headerItem().setText(1, _translate("stdbMainwindow", "dir"))
        self.properties_treeWidget.headerItem().setText(0, _translate("stdbMainwindow", "key"))
        self.properties_treeWidget.headerItem().setText(1, _translate("stdbMainwindow", "value"))
        __sortingEnabled = self.properties_treeWidget.isSortingEnabled()
        self.properties_treeWidget.setSortingEnabled(False)
        self.properties_treeWidget.topLevelItem(0).setText(0, _translate("stdbMainwindow", "Space Group"))
        self.properties_treeWidget.topLevelItem(1).setText(0, _translate("stdbMainwindow", "Unit Cell"))
        self.properties_treeWidget.topLevelItem(1).child(0).setText(0, _translate("stdbMainwindow", "a"))
        self.properties_treeWidget.topLevelItem(1).child(1).setText(0, _translate("stdbMainwindow", "b"))
        self.properties_treeWidget.topLevelItem(1).child(2).setText(0, _translate("stdbMainwindow", "c"))
        self.properties_treeWidget.topLevelItem(1).child(3).setText(0, _translate("stdbMainwindow", "alpha"))
        self.properties_treeWidget.topLevelItem(1).child(4).setText(0, _translate("stdbMainwindow", "beta"))
        self.properties_treeWidget.topLevelItem(1).child(5).setText(0, _translate("stdbMainwindow", "gamma"))
        self.properties_treeWidget.topLevelItem(2).setText(0, _translate("stdbMainwindow", "Temperature"))
        self.properties_treeWidget.topLevelItem(3).setText(0, _translate("stdbMainwindow", "Sum Formula"))
        self.properties_treeWidget.topLevelItem(4).setText(0, _translate("stdbMainwindow", "Z"))
        self.properties_treeWidget.topLevelItem(5).setText(0, _translate("stdbMainwindow", "R1"))
        self.properties_treeWidget.topLevelItem(5).child(0).setText(0, _translate("stdbMainwindow", "I>2s"))
        self.properties_treeWidget.topLevelItem(5).child(1).setText(0, _translate("stdbMainwindow", "All Refl."))
        self.properties_treeWidget.topLevelItem(6).setText(0, _translate("stdbMainwindow", "wR2"))
        self.properties_treeWidget.topLevelItem(6).child(0).setText(0, _translate("stdbMainwindow", "All Refl."))
        self.properties_treeWidget.topLevelItem(6).child(1).setText(0, _translate("stdbMainwindow", "I>2s"))
        self.properties_treeWidget.topLevelItem(7).setText(0, _translate("stdbMainwindow", "Goof"))
        self.properties_treeWidget.topLevelItem(8).setText(0, _translate("stdbMainwindow", "Max. Shift"))
        self.properties_treeWidget.topLevelItem(9).setText(0, _translate("stdbMainwindow", "Peak"))
        self.properties_treeWidget.topLevelItem(10).setText(0, _translate("stdbMainwindow", "Hole"))
        self.properties_treeWidget.setSortingEnabled(__sortingEnabled)
        self.fileMenu.setTitle(_translate("stdbMainwindow", "File"))
        self.menuEdit.setTitle(_translate("stdbMainwindow", "edit"))
        self.menuHelp.setTitle(_translate("stdbMainwindow", "help"))
        self.actionExit.setText(_translate("stdbMainwindow", "Exit"))
        self.actionExit.setShortcut(_translate("stdbMainwindow", "Ctrl+Q"))
        self.actionEdit_dataset.setText(_translate("stdbMainwindow", "edit dataset"))
        self.actionUser_manual.setText(_translate("stdbMainwindow", "user manual"))
        self.actionImport_file.setText(_translate("stdbMainwindow", "Import file"))
        self.actionImport_directory.setText(_translate("stdbMainwindow", "Import directory"))
        self.actionExport_Database_s.setText(_translate("stdbMainwindow", "Export Database(s)"))
        self.actionOptions.setText(_translate("stdbMainwindow", "Options"))

