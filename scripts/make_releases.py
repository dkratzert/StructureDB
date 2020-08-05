#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This script has to be run from the main dir e.g. D:\GitHub\StructureFinder
"""
import os
import sys
from pathlib import Path

app_path = str(Path(os.path.dirname(os.path.abspath(__file__))).parent)
main_path = str(Path(os.path.dirname(os.path.abspath(__file__))))
sys.path.extend([app_path, main_path])
import subprocess

from PyQt5 import uic

from misc.version import VERSION
from scripts.create_zipfile import make_zip, files
from scripts.version_numbers import process_iss, disable_debug, isspath, pypath

print("Updating version numbers to version {} ...".format(VERSION))

# Update version numbers in .iss files
for i in isspath:
    process_iss(i)

# disable all debug variables:
for i in pypath:
    disable_debug(i)

print("Version numbers updated.")

try:
    print(os.path.abspath('./gui'))
    uic.compileUiDir('./gui')
    print('recompiled ui')
except:
    print("Unable to compile UI!")
    raise


def make_distribs():
    # create binary distribution of 64bit variant:
    subprocess.run(['venv/Scripts/pyinstaller.exe',
                    '--clean',
                    '-p',
                    r'D:\Programme\Windows Kits\10\Redist\ucrt\DLLs\x64',
                    r'--add-data=gui;gui',
                    r'--add-data=displaymol;displaymol',
                    r'--add-data=icons;icons',
                    '--hidden-import', 'PyQt5.sip',
                    '-n', 'StructureFinder',
                    '-y', '-i', 'icons/strf.ico',
                    '--windowed',
                    'strf.py'])

    innosetup_compiler = r'C:/Program Files (x86)/Inno Setup 6/ISCC.exe'
    # Run 64bit setup compiler
    subprocess.run([innosetup_compiler,
                    r'scripts\strf-install_win64.iss', ])

    ## Run 32bit setup compiler
    # subprocess.run([innosetup_compiler,
    #                r'scripts\strf-install_win32b.iss', ])


# Make binary distributions:
make_distribs()

# Make a zip file for web interface distribution:
make_zip(files)
