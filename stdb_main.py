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
        stdbMainwindow.resize(958, 722)
        self.centralwidget = QtWidgets.QWidget(stdbMainwindow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.groupBox = QtWidgets.QGroupBox(self.centralwidget)
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
        self.horizontalLayout.addWidget(self.groupBox)
        self.stackedWidget = QtWidgets.QStackedWidget(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(3)
        sizePolicy.setHeightForWidth(self.stackedWidget.sizePolicy().hasHeightForWidth())
        self.stackedWidget.setSizePolicy(sizePolicy)
        self.stackedWidget.setObjectName("stackedWidget")
        self.importResults = QtWidgets.QWidget()
        self.importResults.setObjectName("importResults")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.importResults)
        self.gridLayout_6.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.resultsGroupBox = QtWidgets.QGroupBox(self.importResults)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(4)
        sizePolicy.setHeightForWidth(self.resultsGroupBox.sizePolicy().hasHeightForWidth())
        self.resultsGroupBox.setSizePolicy(sizePolicy)
        self.resultsGroupBox.setObjectName("resultsGroupBox")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.resultsGroupBox)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.formLayout_2 = QtWidgets.QFormLayout()
        self.formLayout_2.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.formLayout_2.setLabelAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.formLayout_2.setObjectName("formLayout_2")
        self.SpaceGrouplabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.SpaceGrouplabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.SpaceGrouplabel.setObjectName("SpaceGrouplabel")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.SpaceGrouplabel)
        self.SpaceGroupLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.SpaceGroupLineEdit.setReadOnly(True)
        self.SpaceGroupLineEdit.setObjectName("SpaceGroupLineEdit")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.SpaceGroupLineEdit)
        self.CellLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.CellLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.CellLabel.setObjectName("CellLabel")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.CellLabel)
        self.aLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.aLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.aLabel.setObjectName("aLabel")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.aLabel)
        self.aLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.aLineEdit.setReadOnly(True)
        self.aLineEdit.setObjectName("aLineEdit")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.aLineEdit)
        self.bLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.bLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.bLabel.setObjectName("bLabel")
        self.formLayout_2.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.bLabel)
        self.bLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.bLineEdit.setReadOnly(True)
        self.bLineEdit.setObjectName("bLineEdit")
        self.formLayout_2.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.bLineEdit)
        self.cLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.cLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.cLabel.setObjectName("cLabel")
        self.formLayout_2.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.cLabel)
        self.cLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.cLineEdit.setReadOnly(True)
        self.cLineEdit.setObjectName("cLineEdit")
        self.formLayout_2.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.cLineEdit)
        self.alphaLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.alphaLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.alphaLabel.setObjectName("alphaLabel")
        self.formLayout_2.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.alphaLabel)
        self.alphaLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.alphaLineEdit.setReadOnly(True)
        self.alphaLineEdit.setObjectName("alphaLineEdit")
        self.formLayout_2.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.alphaLineEdit)
        self.betaLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.betaLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.betaLabel.setObjectName("betaLabel")
        self.formLayout_2.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.betaLabel)
        self.betaLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.betaLineEdit.setObjectName("betaLineEdit")
        self.formLayout_2.setWidget(6, QtWidgets.QFormLayout.FieldRole, self.betaLineEdit)
        self.gammaLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.gammaLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.gammaLabel.setObjectName("gammaLabel")
        self.formLayout_2.setWidget(7, QtWidgets.QFormLayout.LabelRole, self.gammaLabel)
        self.gammaLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.gammaLineEdit.setObjectName("gammaLineEdit")
        self.formLayout_2.setWidget(7, QtWidgets.QFormLayout.FieldRole, self.gammaLineEdit)
        self.zLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.zLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.zLabel.setObjectName("zLabel")
        self.formLayout_2.setWidget(8, QtWidgets.QFormLayout.LabelRole, self.zLabel)
        self.zLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.zLineEdit.setObjectName("zLineEdit")
        self.formLayout_2.setWidget(8, QtWidgets.QFormLayout.FieldRole, self.zLineEdit)
        self.label = QtWidgets.QLabel(self.resultsGroupBox)
        self.label.setObjectName("label")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.label)
        self.horizontalLayout_2.addLayout(self.formLayout_2)
        spacerItem4 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem4)
        self.formLayout_3 = QtWidgets.QFormLayout()
        self.formLayout_3.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.formLayout_3.setLabelAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.formLayout_3.setObjectName("formLayout_3")
        self.temperatureLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.temperatureLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.temperatureLabel.setObjectName("temperatureLabel")
        self.formLayout_3.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.temperatureLabel)
        self.temperatureLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.temperatureLineEdit.setObjectName("temperatureLineEdit")
        self.formLayout_3.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.temperatureLineEdit)
        self.r1Label = QtWidgets.QLabel(self.resultsGroupBox)
        font = QtGui.QFont()
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        self.r1Label.setFont(font)
        self.r1Label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.r1Label.setObjectName("r1Label")
        self.formLayout_3.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.r1Label)
        self.r1LineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.r1LineEdit.setObjectName("r1LineEdit")
        self.formLayout_3.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.r1LineEdit)
        self.goofLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.goofLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.goofLabel.setObjectName("goofLabel")
        self.formLayout_3.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.goofLabel)
        self.goofLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.goofLineEdit.setObjectName("goofLineEdit")
        self.formLayout_3.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.goofLineEdit)
        self.maxShiftLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.maxShiftLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.maxShiftLabel.setObjectName("maxShiftLabel")
        self.formLayout_3.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.maxShiftLabel)
        self.maxShiftLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.maxShiftLineEdit.setObjectName("maxShiftLineEdit")
        self.formLayout_3.setWidget(6, QtWidgets.QFormLayout.FieldRole, self.maxShiftLineEdit)
        self.peakLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.peakLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.peakLabel.setObjectName("peakLabel")
        self.formLayout_3.setWidget(7, QtWidgets.QFormLayout.LabelRole, self.peakLabel)
        self.peakLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.peakLineEdit.setObjectName("peakLineEdit")
        self.formLayout_3.setWidget(7, QtWidgets.QFormLayout.FieldRole, self.peakLineEdit)
        self.rintLabel = QtWidgets.QLabel(self.resultsGroupBox)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.rintLabel.setFont(font)
        self.rintLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.rintLabel.setObjectName("rintLabel")
        self.formLayout_3.setWidget(8, QtWidgets.QFormLayout.LabelRole, self.rintLabel)
        self.rintLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.rintLineEdit.setObjectName("rintLineEdit")
        self.formLayout_3.setWidget(8, QtWidgets.QFormLayout.FieldRole, self.rintLineEdit)
        self.rsigmaLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.rsigmaLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.rsigmaLabel.setObjectName("rsigmaLabel")
        self.formLayout_3.setWidget(9, QtWidgets.QFormLayout.LabelRole, self.rsigmaLabel)
        self.rsigmaLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.rsigmaLineEdit.setObjectName("rsigmaLineEdit")
        self.formLayout_3.setWidget(9, QtWidgets.QFormLayout.FieldRole, self.rsigmaLineEdit)
        self.sumFormulaLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.sumFormulaLineEdit.setObjectName("sumFormulaLineEdit")
        self.formLayout_3.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.sumFormulaLineEdit)
        self.wR2LineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.wR2LineEdit.setObjectName("wR2LineEdit")
        self.formLayout_3.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.wR2LineEdit)
        self.sumFormulaLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.sumFormulaLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.sumFormulaLabel.setObjectName("sumFormulaLabel")
        self.formLayout_3.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.sumFormulaLabel)
        self.wR2Label = QtWidgets.QLabel(self.resultsGroupBox)
        font = QtGui.QFont()
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        self.wR2Label.setFont(font)
        self.wR2Label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.wR2Label.setObjectName("wR2Label")
        self.formLayout_3.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.wR2Label)
        self.horizontalLayout_2.addLayout(self.formLayout_3)
        spacerItem5 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem5)
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setLabelAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.formLayout.setObjectName("formLayout")
        self.reflTotalLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.reflTotalLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.reflTotalLabel.setObjectName("reflTotalLabel")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.reflTotalLabel)
        self.reflTotalLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.reflTotalLineEdit.setObjectName("reflTotalLineEdit")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.reflTotalLineEdit)
        self.numParametersLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.numParametersLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.numParametersLabel.setObjectName("numParametersLabel")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.numParametersLabel)
        self.numParametersLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.numParametersLineEdit.setObjectName("numParametersLineEdit")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.numParametersLineEdit)
        self.dataReflnsLabel = QtWidgets.QLabel(self.resultsGroupBox)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.dataReflnsLabel.setFont(font)
        self.dataReflnsLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.dataReflnsLabel.setObjectName("dataReflnsLabel")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.dataReflnsLabel)
        self.dataReflnsLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.dataReflnsLineEdit.setObjectName("dataReflnsLineEdit")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.dataReflnsLineEdit)
        self.numRestraintsLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.numRestraintsLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.numRestraintsLabel.setObjectName("numRestraintsLabel")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.numRestraintsLabel)
        self.numRestraintsLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.numRestraintsLineEdit.setObjectName("numRestraintsLineEdit")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.numRestraintsLineEdit)
        self.thetaMaxLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.thetaMaxLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.thetaMaxLabel.setObjectName("thetaMaxLabel")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.thetaMaxLabel)
        self.thetaMaxLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.thetaMaxLineEdit.setObjectName("thetaMaxLineEdit")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.thetaMaxLineEdit)
        self.thetaFullLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.thetaFullLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.thetaFullLabel.setObjectName("thetaFullLabel")
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.thetaFullLabel)
        self.thetaFullLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.thetaFullLineEdit.setObjectName("thetaFullLineEdit")
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.thetaFullLineEdit)
        self.dLabel = QtWidgets.QLabel(self.resultsGroupBox)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.dLabel.setFont(font)
        self.dLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.dLabel.setObjectName("dLabel")
        self.formLayout.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.dLabel)
        self.dLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.dLineEdit.setObjectName("dLineEdit")
        self.formLayout.setWidget(6, QtWidgets.QFormLayout.FieldRole, self.dLineEdit)
        self.completeLabel = QtWidgets.QLabel(self.resultsGroupBox)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.completeLabel.setFont(font)
        self.completeLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.completeLabel.setObjectName("completeLabel")
        self.formLayout.setWidget(7, QtWidgets.QFormLayout.LabelRole, self.completeLabel)
        self.completeLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.completeLineEdit.setObjectName("completeLineEdit")
        self.formLayout.setWidget(7, QtWidgets.QFormLayout.FieldRole, self.completeLineEdit)
        self.wavelengthLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.wavelengthLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.wavelengthLabel.setObjectName("wavelengthLabel")
        self.formLayout.setWidget(8, QtWidgets.QFormLayout.LabelRole, self.wavelengthLabel)
        self.wavelengthLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.wavelengthLineEdit.setObjectName("wavelengthLineEdit")
        self.formLayout.setWidget(8, QtWidgets.QFormLayout.FieldRole, self.wavelengthLineEdit)
        self.horizontalLayout_2.addLayout(self.formLayout)
        spacerItem6 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem6)
        self.openglVlayout = QtWidgets.QHBoxLayout()
        self.openglVlayout.setObjectName("openglVlayout")
        self.horizontalLayout_2.addLayout(self.openglVlayout)
        self.gridLayout_6.addWidget(self.resultsGroupBox, 3, 0, 1, 1, QtCore.Qt.AlignHCenter)
        self.cifList_treeWidget = QtWidgets.QTreeWidget(self.importResults)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(3)
        sizePolicy.setHeightForWidth(self.cifList_treeWidget.sizePolicy().hasHeightForWidth())
        self.cifList_treeWidget.setSizePolicy(sizePolicy)
        self.cifList_treeWidget.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.cifList_treeWidget.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.cifList_treeWidget.setObjectName("cifList_treeWidget")
        self.gridLayout_6.addWidget(self.cifList_treeWidget, 2, 0, 1, 1)
        self.searchLayout = QtWidgets.QHBoxLayout()
        self.searchLayout.setSpacing(7)
        self.searchLayout.setObjectName("searchLayout")
        self.searchLabel = QtWidgets.QLabel(self.importResults)
        self.searchLabel.setObjectName("searchLabel")
        self.searchLayout.addWidget(self.searchLabel)
        self.searchLineEDit = QtWidgets.QLineEdit(self.importResults)
        self.searchLineEDit.setObjectName("searchLineEDit")
        self.searchLayout.addWidget(self.searchLineEDit)
        self.gridLayout_6.addLayout(self.searchLayout, 1, 0, 1, 1)
        self.stackedWidget.addWidget(self.importResults)
        self.refinement_edit_page = QtWidgets.QWidget()
        self.refinement_edit_page.setObjectName("refinement_edit_page")
        self.stackedWidget.addWidget(self.refinement_edit_page)
        self.stackedWidgetPage2 = QtWidgets.QWidget()
        self.stackedWidgetPage2.setEnabled(False)
        self.stackedWidgetPage2.setObjectName("stackedWidgetPage2")
        self.stackedWidget.addWidget(self.stackedWidgetPage2)
        self.horizontalLayout.addWidget(self.stackedWidget)
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
        self.resultsGroupBox.setTitle(_translate("stdbMainwindow", "Properties"))
        self.SpaceGrouplabel.setText(_translate("stdbMainwindow", "Space Group"))
        self.CellLabel.setText(_translate("stdbMainwindow", "Unit Cell:"))
        self.aLabel.setText(_translate("stdbMainwindow", "a"))
        self.bLabel.setText(_translate("stdbMainwindow", "b"))
        self.cLabel.setText(_translate("stdbMainwindow", "c"))
        self.alphaLabel.setText(_translate("stdbMainwindow", "alpha"))
        self.betaLabel.setText(_translate("stdbMainwindow", "beta"))
        self.gammaLabel.setText(_translate("stdbMainwindow", "gamma"))
        self.zLabel.setText(_translate("stdbMainwindow", "Z"))
        self.label.setText(_translate("stdbMainwindow", "[Å], deg."))
        self.temperatureLabel.setText(_translate("stdbMainwindow", "Temperature [K]"))
        self.r1Label.setText(_translate("stdbMainwindow", "R1"))
        self.goofLabel.setText(_translate("stdbMainwindow", "Goof"))
        self.maxShiftLabel.setText(_translate("stdbMainwindow", "Max Shift/esd"))
        self.peakLabel.setText(_translate("stdbMainwindow", "Peak / Hole [e/Å^3]"))
        self.rintLabel.setText(_translate("stdbMainwindow", "R(int)"))
        self.rsigmaLabel.setText(_translate("stdbMainwindow", "R(σ)"))
        self.sumFormulaLabel.setText(_translate("stdbMainwindow", "Sum Formula"))
        self.wR2Label.setText(_translate("stdbMainwindow", "wR2"))
        self.reflTotalLabel.setText(_translate("stdbMainwindow", "Total num. Refl."))
        self.numParametersLabel.setText(_translate("stdbMainwindow", "Parameters"))
        self.dataReflnsLabel.setText(_translate("stdbMainwindow", "data/param"))
        self.numRestraintsLabel.setText(_translate("stdbMainwindow", "Restraints"))
        self.thetaMaxLabel.setText(_translate("stdbMainwindow", "θ(max) [°]"))
        self.thetaFullLabel.setText(_translate("stdbMainwindow", "θ(full) [°]"))
        self.dLabel.setText(_translate("stdbMainwindow", "d [Å]"))
        self.completeLabel.setText(_translate("stdbMainwindow", "complete [%]"))
        self.wavelengthLabel.setText(_translate("stdbMainwindow", "Wavelength [Å]"))
        self.cifList_treeWidget.setSortingEnabled(True)
        self.cifList_treeWidget.headerItem().setText(0, _translate("stdbMainwindow", "file"))
        self.cifList_treeWidget.headerItem().setText(1, _translate("stdbMainwindow", "dir"))
        self.cifList_treeWidget.headerItem().setText(2, _translate("stdbMainwindow", "Id"))
        self.searchLabel.setText(_translate("stdbMainwindow", "Cell Search:"))
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

