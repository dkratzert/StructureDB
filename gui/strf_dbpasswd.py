# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'D:\GitHub\StructureFinder\./gui\strf_dbpasswd.ui'
#
# Created: Wed May 29 14:10:51 2019
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_PasswdDialog(object):
    def setupUi(self, PasswdDialog):
        PasswdDialog.setObjectName(_fromUtf8("PasswdDialog"))
        PasswdDialog.setWindowModality(QtCore.Qt.ApplicationModal)
        PasswdDialog.resize(400, 213)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(PasswdDialog.sizePolicy().hasHeightForWidth())
        PasswdDialog.setSizePolicy(sizePolicy)
        PasswdDialog.setSizeGripEnabled(False)
        PasswdDialog.setModal(True)
        self.gridLayout = QtGui.QGridLayout(PasswdDialog)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.formLayout = QtGui.QFormLayout()
        self.formLayout.setSizeConstraint(QtGui.QLayout.SetDefaultConstraint)
        self.formLayout.setFieldGrowthPolicy(QtGui.QFormLayout.ExpandingFieldsGrow)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.WarnLabel = QtGui.QLabel(PasswdDialog)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.WarnLabel.setFont(font)
        self.WarnLabel.setObjectName(_fromUtf8("WarnLabel"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.SpanningRole, self.WarnLabel)
        self.Asklabel = QtGui.QLabel(PasswdDialog)
        self.Asklabel.setObjectName(_fromUtf8("Asklabel"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.SpanningRole, self.Asklabel)
        self.NameLabel = QtGui.QLabel(PasswdDialog)
        self.NameLabel.setObjectName(_fromUtf8("NameLabel"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.LabelRole, self.NameLabel)
        self.userNameLineEdit = QtGui.QLineEdit(PasswdDialog)
        self.userNameLineEdit.setObjectName(_fromUtf8("userNameLineEdit"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.FieldRole, self.userNameLineEdit)
        self.PasswordLabel = QtGui.QLabel(PasswdDialog)
        self.PasswordLabel.setObjectName(_fromUtf8("PasswordLabel"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.LabelRole, self.PasswordLabel)
        self.PasswordLineEdit = QtGui.QLineEdit(PasswdDialog)
        self.PasswordLineEdit.setObjectName(_fromUtf8("PasswordLineEdit"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.FieldRole, self.PasswordLineEdit)
        self.IPLabel = QtGui.QLabel(PasswdDialog)
        self.IPLabel.setObjectName(_fromUtf8("IPLabel"))
        self.formLayout.setWidget(4, QtGui.QFormLayout.LabelRole, self.IPLabel)
        self.IPlineEdit = QtGui.QLineEdit(PasswdDialog)
        self.IPlineEdit.setObjectName(_fromUtf8("IPlineEdit"))
        self.formLayout.setWidget(4, QtGui.QFormLayout.FieldRole, self.IPlineEdit)
        self.verticalLayout.addLayout(self.formLayout)
        self.buttonBox = QtGui.QDialogButtonBox(PasswdDialog)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.verticalLayout.addWidget(self.buttonBox)
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Preferred)
        self.verticalLayout.addItem(spacerItem)
        self.gridLayout.addLayout(self.verticalLayout, 0, 0, 1, 1)

        self.retranslateUi(PasswdDialog)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), PasswdDialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), PasswdDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(PasswdDialog)

    def retranslateUi(self, PasswdDialog):
        PasswdDialog.setWindowTitle(_translate("PasswdDialog", "StructureFinder", None))
        self.WarnLabel.setText(_translate("PasswdDialog", "Unable to connect to the APEX database. ", None))
        self.Asklabel.setText(_translate("PasswdDialog", "Please give the correct IP Adress, Username and Password.", None))
        self.NameLabel.setText(_translate("PasswdDialog", "Username", None))
        self.userNameLineEdit.setText(_translate("PasswdDialog", "BrukerPostgreSQL", None))
        self.PasswordLabel.setText(_translate("PasswdDialog", "Password", None))
        self.PasswordLineEdit.setText(_translate("PasswdDialog", "Bruker-PostgreSQL", None))
        self.IPLabel.setText(_translate("PasswdDialog", "IP Adress", None))
        self.IPlineEdit.setText(_translate("PasswdDialog", "localhost", None))

