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
        stdbMainwindow.resize(1064, 781)
        self.centralwidget = QtWidgets.QWidget(stdbMainwindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout_2.setObjectName("gridLayout_2")
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
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        self.gridLayout.addItem(spacerItem, 7, 0, 1, 1)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem1, 9, 0, 1, 1)
        spacerItem2 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        self.gridLayout.addItem(spacerItem2, 5, 0, 1, 1)
        self.importDirButton = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.importDirButton.sizePolicy().hasHeightForWidth())
        self.importDirButton.setSizePolicy(sizePolicy)
        self.importDirButton.setObjectName("importDirButton")
        self.gridLayout.addWidget(self.importDirButton, 0, 0, 1, 1)
        self.openApexDBButton = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.openApexDBButton.sizePolicy().hasHeightForWidth())
        self.openApexDBButton.setSizePolicy(sizePolicy)
        self.openApexDBButton.setObjectName("openApexDBButton")
        self.gridLayout.addWidget(self.openApexDBButton, 6, 0, 1, 1)
        self.importDatabaseButton = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.importDatabaseButton.sizePolicy().hasHeightForWidth())
        self.importDatabaseButton.setSizePolicy(sizePolicy)
        self.importDatabaseButton.setObjectName("importDatabaseButton")
        self.gridLayout.addWidget(self.importDatabaseButton, 2, 0, 1, 1)
        spacerItem3 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        self.gridLayout.addItem(spacerItem3, 1, 0, 1, 1)
        self.saveDatabaseButton = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.saveDatabaseButton.sizePolicy().hasHeightForWidth())
        self.saveDatabaseButton.setSizePolicy(sizePolicy)
        self.saveDatabaseButton.setObjectName("saveDatabaseButton")
        self.gridLayout.addWidget(self.saveDatabaseButton, 8, 0, 1, 1)
        self.closeDatabaseButton = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.closeDatabaseButton.sizePolicy().hasHeightForWidth())
        self.closeDatabaseButton.setSizePolicy(sizePolicy)
        self.closeDatabaseButton.setObjectName("closeDatabaseButton")
        self.gridLayout.addWidget(self.closeDatabaseButton, 10, 0, 1, 1)
        self.gridLayout_2.addWidget(self.groupBox, 0, 0, 1, 1)
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setAccessibleName("")
        self.tabWidget.setAccessibleDescription("")
        self.tabWidget.setAutoFillBackground(False)
        self.tabWidget.setStyleSheet("")
        self.tabWidget.setObjectName("tabWidget")
        self.maintab = QtWidgets.QWidget()
        self.maintab.setObjectName("maintab")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.maintab)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.importResults = QtWidgets.QGridLayout()
        self.importResults.setSizeConstraint(QtWidgets.QLayout.SetMaximumSize)
        self.importResults.setObjectName("importResults")
        self.cifList_treeWidget = QtWidgets.QTreeWidget(self.maintab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(3)
        sizePolicy.setHeightForWidth(self.cifList_treeWidget.sizePolicy().hasHeightForWidth())
        self.cifList_treeWidget.setSizePolicy(sizePolicy)
        self.cifList_treeWidget.setAutoFillBackground(False)
        self.cifList_treeWidget.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.cifList_treeWidget.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.cifList_treeWidget.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.cifList_treeWidget.setObjectName("cifList_treeWidget")
        self.importResults.addWidget(self.cifList_treeWidget, 1, 0, 1, 1)
        self.searchLayout = QtWidgets.QHBoxLayout()
        self.searchLayout.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.searchLayout.setSpacing(7)
        self.searchLayout.setObjectName("searchLayout")
        self.searchLabel = QtWidgets.QLabel(self.maintab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.searchLabel.sizePolicy().hasHeightForWidth())
        self.searchLabel.setSizePolicy(sizePolicy)
        self.searchLabel.setObjectName("searchLabel")
        self.searchLayout.addWidget(self.searchLabel)
        self.searchCellLineEDit = QtWidgets.QLineEdit(self.maintab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(3)
        sizePolicy.setHeightForWidth(self.searchCellLineEDit.sizePolicy().hasHeightForWidth())
        self.searchCellLineEDit.setSizePolicy(sizePolicy)
        self.searchCellLineEDit.setObjectName("searchCellLineEDit")
        self.searchLayout.addWidget(self.searchCellLineEDit)
        self.moreResultsCheckBox = QtWidgets.QCheckBox(self.maintab)
        self.moreResultsCheckBox.setObjectName("moreResultsCheckBox")
        self.searchLayout.addWidget(self.moreResultsCheckBox)
        spacerItem4 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        self.searchLayout.addItem(spacerItem4)
        self.line = QtWidgets.QFrame(self.maintab)
        self.line.setFrameShape(QtWidgets.QFrame.VLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.searchLayout.addWidget(self.line)
        self.txtSearchLabel = QtWidgets.QLabel(self.maintab)
        self.txtSearchLabel.setObjectName("txtSearchLabel")
        self.searchLayout.addWidget(self.txtSearchLabel)
        self.txtSearchEdit = QtWidgets.QLineEdit(self.maintab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(2)
        sizePolicy.setHeightForWidth(self.txtSearchEdit.sizePolicy().hasHeightForWidth())
        self.txtSearchEdit.setSizePolicy(sizePolicy)
        self.txtSearchEdit.setObjectName("txtSearchEdit")
        self.searchLayout.addWidget(self.txtSearchEdit)
        self.importResults.addLayout(self.searchLayout, 0, 0, 1, 1)
        spacerItem5 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.importResults.addItem(spacerItem5, 3, 0, 1, 1)
        self.resultsGroupBox = QtWidgets.QGroupBox(self.maintab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(4)
        sizePolicy.setVerticalStretch(4)
        sizePolicy.setHeightForWidth(self.resultsGroupBox.sizePolicy().hasHeightForWidth())
        self.resultsGroupBox.setSizePolicy(sizePolicy)
        self.resultsGroupBox.setStyleSheet("")
        self.resultsGroupBox.setFlat(False)
        self.resultsGroupBox.setObjectName("resultsGroupBox")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.resultsGroupBox)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout()
        self.verticalLayout_6.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.groupBox_4 = QtWidgets.QGroupBox(self.resultsGroupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_4.sizePolicy().hasHeightForWidth())
        self.groupBox_4.setSizePolicy(sizePolicy)
        self.groupBox_4.setMinimumSize(QtCore.QSize(210, 0))
        self.groupBox_4.setAutoFillBackground(False)
        self.groupBox_4.setTitle("")
        self.groupBox_4.setFlat(False)
        self.groupBox_4.setObjectName("groupBox_4")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.groupBox_4)
        self.gridLayout_5.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_5.setHorizontalSpacing(6)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.cellField = QtWidgets.QLabel(self.groupBox_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.cellField.sizePolicy().hasHeightForWidth())
        self.cellField.setSizePolicy(sizePolicy)
        self.cellField.setMinimumSize(QtCore.QSize(0, 75))
        self.cellField.setBaseSize(QtCore.QSize(160, 75))
        self.cellField.setAutoFillBackground(False)
        self.cellField.setStyleSheet("border-color: rgb(53, 53, 53);")
        self.cellField.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.cellField.setFrameShadow(QtWidgets.QFrame.Plain)
        self.cellField.setText("")
        self.cellField.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.cellField.setObjectName("cellField")
        self.gridLayout_5.addWidget(self.cellField, 1, 1, 1, 1)
        self.CellLabel = QtWidgets.QLabel(self.groupBox_4)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.CellLabel.setFont(font)
        self.CellLabel.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.CellLabel.setObjectName("CellLabel")
        self.gridLayout_5.addWidget(self.CellLabel, 0, 1, 1, 1)
        self.verticalLayout_6.addWidget(self.groupBox_4)
        self.label = QtWidgets.QLabel(self.resultsGroupBox)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.verticalLayout_6.addWidget(self.label)
        self.formLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.formLabel.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.formLabel.sizePolicy().hasHeightForWidth())
        self.formLabel.setSizePolicy(sizePolicy)
        self.formLabel.setMinimumSize(QtCore.QSize(0, 21))
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setKerning(False)
        self.formLabel.setFont(font)
        self.formLabel.setAutoFillBackground(False)
        self.formLabel.setStyleSheet("background-color:  rgba(255, 255, 255, 255);")
        self.formLabel.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.formLabel.setLineWidth(1)
        self.formLabel.setText("")
        self.formLabel.setTextFormat(QtCore.Qt.RichText)
        self.formLabel.setWordWrap(True)
        self.formLabel.setObjectName("formLabel")
        self.verticalLayout_6.addWidget(self.formLabel)
        spacerItem6 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_6.addItem(spacerItem6)
        spacerItem7 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_6.addItem(spacerItem7)
        self.horizontalLayout_2.addLayout(self.verticalLayout_6)
        spacerItem8 = QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem8)
        self.formLayout_3 = QtWidgets.QFormLayout()
        self.formLayout_3.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)
        self.formLayout_3.setLabelAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.formLayout_3.setObjectName("formLayout_3")
        self.SpaceGrouplabel = QtWidgets.QLabel(self.resultsGroupBox)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.SpaceGrouplabel.setFont(font)
        self.SpaceGrouplabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.SpaceGrouplabel.setObjectName("SpaceGrouplabel")
        self.formLayout_3.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.SpaceGrouplabel)
        self.SpaceGroupLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.SpaceGroupLineEdit.setReadOnly(True)
        self.SpaceGroupLineEdit.setObjectName("SpaceGroupLineEdit")
        self.formLayout_3.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.SpaceGroupLineEdit)
        self.zLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.zLabel.setObjectName("zLabel")
        self.formLayout_3.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.zLabel)
        self.zLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.zLineEdit.setReadOnly(True)
        self.zLineEdit.setObjectName("zLineEdit")
        self.formLayout_3.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.zLineEdit)
        self.temperatureLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.temperatureLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.temperatureLabel.setObjectName("temperatureLabel")
        self.formLayout_3.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.temperatureLabel)
        self.temperatureLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.temperatureLineEdit.setReadOnly(True)
        self.temperatureLineEdit.setObjectName("temperatureLineEdit")
        self.formLayout_3.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.temperatureLineEdit)
        self.wR2Label = QtWidgets.QLabel(self.resultsGroupBox)
        font = QtGui.QFont()
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        self.wR2Label.setFont(font)
        self.wR2Label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.wR2Label.setObjectName("wR2Label")
        self.formLayout_3.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.wR2Label)
        self.wR2LineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.wR2LineEdit.setReadOnly(True)
        self.wR2LineEdit.setObjectName("wR2LineEdit")
        self.formLayout_3.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.wR2LineEdit)
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
        self.r1LineEdit.setReadOnly(True)
        self.r1LineEdit.setObjectName("r1LineEdit")
        self.formLayout_3.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.r1LineEdit)
        self.goofLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.goofLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.goofLabel.setObjectName("goofLabel")
        self.formLayout_3.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.goofLabel)
        self.goofLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.goofLineEdit.setReadOnly(True)
        self.goofLineEdit.setObjectName("goofLineEdit")
        self.formLayout_3.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.goofLineEdit)
        self.maxShiftLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.maxShiftLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.maxShiftLabel.setObjectName("maxShiftLabel")
        self.formLayout_3.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.maxShiftLabel)
        self.maxShiftLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.maxShiftLineEdit.setReadOnly(True)
        self.maxShiftLineEdit.setObjectName("maxShiftLineEdit")
        self.formLayout_3.setWidget(6, QtWidgets.QFormLayout.FieldRole, self.maxShiftLineEdit)
        self.peakLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.peakLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.peakLabel.setObjectName("peakLabel")
        self.formLayout_3.setWidget(7, QtWidgets.QFormLayout.LabelRole, self.peakLabel)
        self.peakLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.peakLineEdit.setReadOnly(True)
        self.peakLineEdit.setObjectName("peakLineEdit")
        self.formLayout_3.setWidget(7, QtWidgets.QFormLayout.FieldRole, self.peakLineEdit)
        self.rintLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.rintLineEdit.setReadOnly(True)
        self.rintLineEdit.setObjectName("rintLineEdit")
        self.formLayout_3.setWidget(8, QtWidgets.QFormLayout.FieldRole, self.rintLineEdit)
        self.rsigmaLabel = QtWidgets.QLabel(self.resultsGroupBox)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.rsigmaLabel.setFont(font)
        self.rsigmaLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.rsigmaLabel.setObjectName("rsigmaLabel")
        self.formLayout_3.setWidget(9, QtWidgets.QFormLayout.LabelRole, self.rsigmaLabel)
        self.rsigmaLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.rsigmaLineEdit.setReadOnly(True)
        self.rsigmaLineEdit.setObjectName("rsigmaLineEdit")
        self.formLayout_3.setWidget(9, QtWidgets.QFormLayout.FieldRole, self.rsigmaLineEdit)
        self.rintLabel = QtWidgets.QLabel(self.resultsGroupBox)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.rintLabel.setFont(font)
        self.rintLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.rintLabel.setObjectName("rintLabel")
        self.formLayout_3.setWidget(8, QtWidgets.QFormLayout.LabelRole, self.rintLabel)
        self.horizontalLayout_2.addLayout(self.formLayout_3)
        spacerItem9 = QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem9)
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setLabelAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.formLayout.setObjectName("formLayout")
        self.reflTotalLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.reflTotalLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.reflTotalLabel.setObjectName("reflTotalLabel")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.reflTotalLabel)
        self.reflTotalLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.reflTotalLineEdit.setReadOnly(True)
        self.reflTotalLineEdit.setObjectName("reflTotalLineEdit")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.reflTotalLineEdit)
        self.numParametersLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.numParametersLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.numParametersLabel.setObjectName("numParametersLabel")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.numParametersLabel)
        self.numParametersLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.numParametersLineEdit.setReadOnly(True)
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
        self.dataReflnsLineEdit.setReadOnly(True)
        self.dataReflnsLineEdit.setObjectName("dataReflnsLineEdit")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.dataReflnsLineEdit)
        self.numRestraintsLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.numRestraintsLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.numRestraintsLabel.setObjectName("numRestraintsLabel")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.numRestraintsLabel)
        self.numRestraintsLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.numRestraintsLineEdit.setReadOnly(True)
        self.numRestraintsLineEdit.setObjectName("numRestraintsLineEdit")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.numRestraintsLineEdit)
        self.thetaMaxLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.thetaMaxLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.thetaMaxLabel.setObjectName("thetaMaxLabel")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.thetaMaxLabel)
        self.thetaMaxLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.thetaMaxLineEdit.setReadOnly(True)
        self.thetaMaxLineEdit.setObjectName("thetaMaxLineEdit")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.thetaMaxLineEdit)
        self.thetaFullLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.thetaFullLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.thetaFullLabel.setObjectName("thetaFullLabel")
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.thetaFullLabel)
        self.thetaFullLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.thetaFullLineEdit.setReadOnly(True)
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
        self.dLineEdit.setReadOnly(True)
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
        self.completeLineEdit.setReadOnly(True)
        self.completeLineEdit.setObjectName("completeLineEdit")
        self.formLayout.setWidget(7, QtWidgets.QFormLayout.FieldRole, self.completeLineEdit)
        self.wavelengthLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.wavelengthLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.wavelengthLabel.setObjectName("wavelengthLabel")
        self.formLayout.setWidget(8, QtWidgets.QFormLayout.LabelRole, self.wavelengthLabel)
        self.wavelengthLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.wavelengthLineEdit.setReadOnly(True)
        self.wavelengthLineEdit.setObjectName("wavelengthLineEdit")
        self.formLayout.setWidget(8, QtWidgets.QFormLayout.FieldRole, self.wavelengthLineEdit)
        self.cCDCNumberLabel = QtWidgets.QLabel(self.resultsGroupBox)
        self.cCDCNumberLabel.setObjectName("cCDCNumberLabel")
        self.formLayout.setWidget(9, QtWidgets.QFormLayout.LabelRole, self.cCDCNumberLabel)
        self.cCDCNumberLineEdit = QtWidgets.QLineEdit(self.resultsGroupBox)
        self.cCDCNumberLineEdit.setReadOnly(True)
        self.cCDCNumberLineEdit.setObjectName("cCDCNumberLineEdit")
        self.formLayout.setWidget(9, QtWidgets.QFormLayout.FieldRole, self.cCDCNumberLineEdit)
        self.horizontalLayout_2.addLayout(self.formLayout)
        self.openglview = QtWidgets.QGroupBox(self.resultsGroupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.openglview.sizePolicy().hasHeightForWidth())
        self.openglview.setSizePolicy(sizePolicy)
        self.openglview.setAlignment(QtCore.Qt.AlignCenter)
        self.openglview.setFlat(False)
        self.openglview.setCheckable(False)
        self.openglview.setObjectName("openglview")
        self.ogllayout = QtWidgets.QVBoxLayout(self.openglview)
        self.ogllayout.setContentsMargins(-1, -1, 0, 0)
        self.ogllayout.setSpacing(0)
        self.ogllayout.setObjectName("ogllayout")
        self.horizontalLayout_2.addWidget(self.openglview)
        self.importResults.addWidget(self.resultsGroupBox, 2, 0, 1, 1)
        self.gridLayout_3.addLayout(self.importResults, 0, 0, 1, 1)
        self.tabWidget.addTab(self.maintab, "")
        self.searchtab = QtWidgets.QWidget()
        self.searchtab.setObjectName("searchtab")
        self.tabWidget.addTab(self.searchtab, "")
        self.optiontab = QtWidgets.QWidget()
        self.optiontab.setObjectName("optiontab")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.optiontab)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.groupBox_3 = QtWidgets.QGroupBox(self.optiontab)
        self.groupBox_3.setObjectName("groupBox_3")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.groupBox_3)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout_3.addLayout(self.verticalLayout_2)
        self.horizontalLayout.addWidget(self.groupBox_3)
        self.groupBox_2 = QtWidgets.QGroupBox(self.optiontab)
        self.groupBox_2.setObjectName("groupBox_2")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.groupBox_2)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout()
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.verticalLayout_5.addLayout(self.verticalLayout_4)
        self.horizontalLayout.addWidget(self.groupBox_2)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.gridLayout_4.addLayout(self.horizontalLayout, 0, 0, 1, 1)
        self.tabWidget.addTab(self.optiontab, "")
        self.allEntrysTab = QtWidgets.QWidget()
        self.allEntrysTab.setObjectName("allEntrysTab")
        self.verticalLayout_8 = QtWidgets.QVBoxLayout(self.allEntrysTab)
        self.verticalLayout_8.setObjectName("verticalLayout_8")
        self.verticalLayout_7 = QtWidgets.QVBoxLayout()
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.allCifTreeWidget = QtWidgets.QTreeWidget(self.allEntrysTab)
        self.allCifTreeWidget.setObjectName("allCifTreeWidget")
        self.verticalLayout_7.addWidget(self.allCifTreeWidget)
        self.verticalLayout_8.addLayout(self.verticalLayout_7)
        self.tabWidget.addTab(self.allEntrysTab, "")
        self.gridLayout_2.addWidget(self.tabWidget, 0, 1, 1, 1)
        stdbMainwindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(stdbMainwindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1064, 21))
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
        self.actionSave_Database = QtWidgets.QAction(stdbMainwindow)
        self.actionSave_Database.setObjectName("actionSave_Database")
        self.actionOptions = QtWidgets.QAction(stdbMainwindow)
        self.actionOptions.setObjectName("actionOptions")
        self.actionAdvanced_Search = QtWidgets.QAction(stdbMainwindow)
        self.actionAdvanced_Search.setObjectName("actionAdvanced_Search")
        self.actionAbout_StructureDB = QtWidgets.QAction(stdbMainwindow)
        self.actionAbout_StructureDB.setObjectName("actionAbout_StructureDB")
        self.actionClose_Database = QtWidgets.QAction(stdbMainwindow)
        self.actionClose_Database.setObjectName("actionClose_Database")
        self.fileMenu.addAction(self.actionImport_directory)
        self.fileMenu.addAction(self.actionImport_file)
        self.fileMenu.addAction(self.actionSave_Database)
        self.fileMenu.addAction(self.actionClose_Database)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.actionExit)
        self.menuEdit.addSeparator()
        self.menuEdit.addAction(self.actionOptions)
        self.menuHelp.addAction(self.actionUser_manual)
        self.menuHelp.addAction(self.actionAbout_StructureDB)
        self.menubar.addAction(self.fileMenu.menuAction())
        self.menubar.addAction(self.menuEdit.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(stdbMainwindow)
        self.tabWidget.setCurrentIndex(0)
        self.actionExit.triggered.connect(stdbMainwindow.close)
        QtCore.QMetaObject.connectSlotsByName(stdbMainwindow)

    def retranslateUi(self, stdbMainwindow):
        _translate = QtCore.QCoreApplication.translate
        stdbMainwindow.setWindowTitle(_translate("stdbMainwindow", "StructureFinder"))
        self.importDirButton.setText(_translate("stdbMainwindow", "Import Directory"))
        self.openApexDBButton.setText(_translate("stdbMainwindow", "Open APEX Database"))
        self.importDatabaseButton.setText(_translate("stdbMainwindow", "Open Database File"))
        self.saveDatabaseButton.setText(_translate("stdbMainwindow", "Save Database"))
        self.closeDatabaseButton.setText(_translate("stdbMainwindow", "Close Database"))
        self.cifList_treeWidget.setSortingEnabled(True)
        self.cifList_treeWidget.headerItem().setText(0, _translate("stdbMainwindow", "Filename"))
        self.cifList_treeWidget.headerItem().setText(1, _translate("stdbMainwindow", "Data Name"))
        self.cifList_treeWidget.headerItem().setText(2, _translate("stdbMainwindow", "Directory"))
        self.cifList_treeWidget.headerItem().setText(3, _translate("stdbMainwindow", "Id"))
        self.searchLabel.setText(_translate("stdbMainwindow", "Cell Search:"))
        self.moreResultsCheckBox.setText(_translate("stdbMainwindow", "More results"))
        self.txtSearchLabel.setText(_translate("stdbMainwindow", "Search text:"))
        self.resultsGroupBox.setTitle(_translate("stdbMainwindow", "Properties"))
        self.CellLabel.setText(_translate("stdbMainwindow", "Unit Cell:"))
        self.label.setText(_translate("stdbMainwindow", "Sum Formula:"))
        self.SpaceGrouplabel.setText(_translate("stdbMainwindow", "Space Group"))
        self.zLabel.setText(_translate("stdbMainwindow", "Z"))
        self.temperatureLabel.setText(_translate("stdbMainwindow", "Temperature [K]"))
        self.wR2Label.setText(_translate("stdbMainwindow", "wR2"))
        self.r1Label.setText(_translate("stdbMainwindow", "R1"))
        self.goofLabel.setText(_translate("stdbMainwindow", "Goof"))
        self.maxShiftLabel.setText(_translate("stdbMainwindow", "Max Shift/esd"))
        self.peakLabel.setText(_translate("stdbMainwindow", "<html><head/><body><p>Peak / Hole [eÅ<span style=\" vertical-align:super;\">-3</span>]</p></body></html>"))
        self.rsigmaLabel.setText(_translate("stdbMainwindow", "R(σ)"))
        self.rintLabel.setText(_translate("stdbMainwindow", "R(int)"))
        self.reflTotalLabel.setText(_translate("stdbMainwindow", "Total num. Refl."))
        self.numParametersLabel.setText(_translate("stdbMainwindow", "Parameters"))
        self.dataReflnsLabel.setText(_translate("stdbMainwindow", "data/param"))
        self.numRestraintsLabel.setText(_translate("stdbMainwindow", "Restraints"))
        self.thetaMaxLabel.setText(_translate("stdbMainwindow", "θ(max) [°]"))
        self.thetaFullLabel.setText(_translate("stdbMainwindow", "θ(full) [°]"))
        self.dLabel.setText(_translate("stdbMainwindow", "d [Å]"))
        self.completeLabel.setText(_translate("stdbMainwindow", "complete [%]"))
        self.wavelengthLabel.setText(_translate("stdbMainwindow", "Wavelength [Å]"))
        self.cCDCNumberLabel.setText(_translate("stdbMainwindow", "CCDC Number"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.maintab), _translate("stdbMainwindow", "Structures"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.searchtab), _translate("stdbMainwindow", "Advanced Search"))
        self.groupBox_3.setTitle(_translate("stdbMainwindow", "GroupBox"))
        self.groupBox_2.setTitle(_translate("stdbMainwindow", "GroupBox"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.optiontab), _translate("stdbMainwindow", "Options"))
        self.allCifTreeWidget.headerItem().setText(0, _translate("stdbMainwindow", "Keys"))
        self.allCifTreeWidget.headerItem().setText(1, _translate("stdbMainwindow", "Values"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.allEntrysTab), _translate("stdbMainwindow", "All cif values"))
        self.fileMenu.setTitle(_translate("stdbMainwindow", "File"))
        self.menuEdit.setTitle(_translate("stdbMainwindow", "edit"))
        self.menuHelp.setTitle(_translate("stdbMainwindow", "help"))
        self.actionExit.setText(_translate("stdbMainwindow", "Exit Program"))
        self.actionExit.setShortcut(_translate("stdbMainwindow", "Ctrl+Q"))
        self.actionEdit_dataset.setText(_translate("stdbMainwindow", "edit dataset"))
        self.actionUser_manual.setText(_translate("stdbMainwindow", "User Manual"))
        self.actionImport_file.setText(_translate("stdbMainwindow", "Open Database"))
        self.actionImport_directory.setText(_translate("stdbMainwindow", "Import Directory"))
        self.actionSave_Database.setText(_translate("stdbMainwindow", "Save Database"))
        self.actionOptions.setText(_translate("stdbMainwindow", "Options"))
        self.actionAdvanced_Search.setText(_translate("stdbMainwindow", "Advanced Search"))
        self.actionAbout_StructureDB.setText(_translate("stdbMainwindow", "About StructureDB"))
        self.actionClose_Database.setText(_translate("stdbMainwindow", "Close Database"))

