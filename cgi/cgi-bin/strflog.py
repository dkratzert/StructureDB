#!/usr/local/bin/python3.6
##!C:\tools\Python-3.6.2_64\python.exe

import pathlib
import pprint
from string import Template
from urllib import parse

import cgi
from lattice import lattice
from pymatgen.core import mat_lattice
from searcher import database_handler

import cgitb
cgitb.enable()

from searcher.database_handler import StructureTable


def application():
    """
    The main application of the StructureFinder web interface.
    """
    print("Content-type:text/html\r\n\r\n")
    print('<html>')
    print('<head>')
    print('<title>Hello Word - First CGI Program</title>')
    print('</head>')
    print('<body>')
    print('<h2>Hello Word! This is my first CGI program</h2>')
    print('</body>')
    print('</html>')


def find_cell(structures: StructureTable, cellstr: str) -> list:
    """
    Finds unit cells in db. Rsturns hits a a list of ids.
    """
    try:
        cell = [float(x) for x in cellstr.strip().split()]
    except (TypeError, ValueError) as e:
        print(e)
        return []
    if len(cell) != 6:
        print("No valid cell!")
        return []
    # if self.ui.moreResultsCheckBox.isChecked() or \
    #        self.ui.ad_moreResultscheckBox.isChecked():
    threshold = 0.08
    ltol = 0.09
    atol = 1.8
    # else:
    #    threshold = 0.03
    #    ltol = 0.001
    #    atol = 1
    # try:
    volume = lattice.vol_unitcell(*cell)
    idlist = structures.find_by_volume(volume, threshold)
    idlist2 = []
    if idlist:
        lattice1 = mat_lattice.Lattice.from_parameters_niggli_reduced(*cell)
        for num, i in enumerate(idlist):
            request = """select * from cell where StructureId = {}""".format(i)
            dic = structures.get_row_as_dict(request)
            try:
                lattice2 = mat_lattice.Lattice.from_parameters(
                    float(dic['a']),
                    float(dic['b']),
                    float(dic['c']),
                    float(dic['alpha']),
                    float(dic['beta']),
                    float(dic['gamma']))
            except ValueError:
                continue
            map = lattice1.find_mapping(lattice2, ltol, atol, skip_rotation_matrix=True)
            if map:
                idlist2.append(i)
    if idlist2:
        return idlist2


def search_text(structures: StructureTable, search_string: str) -> list:
    """
    searches db for given text
    """
    idlist = []
    if len(search_string) == 0:
        return []
    if len(search_string) >= 2:
        if "*" not in search_string:
            search_string = "{}{}{}".format('*', search_string, '*')
    try:
        #  bad hack, should make this return ids like cell search
        idlist = [x[0] for x in structures.find_by_strings(search_string)]
    except AttributeError as e:
        print("Error 1")
        print(e)
    return idlist


def process_data(structures: StructureTable, idlist: list=None):
    """
    Structure.Id,           0
    Structure.measurement,  1
    Structure.path,         2
    Structure.filename,     3
    Structure.dataname      4
    """
    print("process data ###")
    if not structures:
        return []
    table_string = ""
    for i in structures.get_all_structure_names(idlist):
        table_string += '<tr id={3}> <td> {0} </a></td> ' \
                        '     <td> {1} </a></td> ' \
                        '     <td> {2} </a></td> </tr> \n' \
                            .format(i[3].decode('utf-8', errors='ignore'),
                                    i[4].decode('utf-8', errors='ignore'),
                                    i[2].decode('utf-8', errors='ignore'),
                                    i[0]
                                    )
        # i[0] -> id
    p = pathlib.Path("./cgi/strflog_Template.htm")
    t = Template(p.read_bytes().decode('utf-8', 'ignore'))
    replacedict = {"logtablecolumns": table_string, "CSearch": "Search", "TSearch": "Search"}
    return str(t.safe_substitute(replacedict)).encode('utf-8', 'ignore')


if __name__ == "__main__":
    application()
    """
    try:
        import wsgiref.simple_server

        server = wsgiref.simple_server.make_server('127.0.0.1', 8000, application)
        server.serve_forever()
        print("Webserver running...")
    except KeyboardInterrupt:
        print("Webserver stopped...")
    """