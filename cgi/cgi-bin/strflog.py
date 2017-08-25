#!C:\tools\Python-3.6.2_64\python.exe
##!/usr/local/bin/python3

import pathlib
from string import Template
from urllib import parse
from wsgiref.util import setup_testing_defaults

from lattice import lattice
from searcher import database_handler


# cgitb.enable(display=1, logdir="./log")


def application(environ, start_response):
    """
    The main application of the StructureFinder web interface.
    """
    setup_testing_defaults(environ)
    # the environment variable CONTENT_LENGTH may be empty or missing
    try:
        request_body_size = int(environ.get('CONTENT_LENGTH', 0))
    except ValueError:
        request_body_size = 0
    # When the method is POST the variable will be sent
    # in the HTTP request body which is passed by the WSGI server
    # in the file like wsgi.input environment variable.
    # pprint.pprint(environ)
    request_body = environ['wsgi.input'].read(request_body_size).decode('utf-8')
    d = parse.parse_qs(request_body)
    if d.get("cell"):
        cell = (d.get("cell"))
    # pprint.pprint('request_body:', d)
    status = '200 OK'
    headers = [('Content-type', 'text/html; charset=utf-8')]
    start_response(status, headers)
    dbfilename = "./structuredb.sqlite"
    structures = database_handler.StructureTable(dbfilename)
    txt = process_data(structures)
    return [txt]


def find_cell(cellstr: str):
    """
    Finds unit cells in db. Rsturns hits a a list of ids.
    """
    try:
        cell = [float(x) for x in cellstr.strip().split()]
    except (TypeError, ValueError):
        return []
    if len(cell) != 6:
        # self.statusBar().showMessage('Not a valid unit cell!', msecs=3000)
        # self.show_full_list()
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
    # idlist = structures.find_by_volume(volume, threshold)


def process_data(structures):
    """
    Structure.Id             0
    Structure.measurement    1
    Structure.path           2
    Structure.filename       3
    Structure.dataname       4
    """
    if not structures:
        return []
    table_string = ""
    for i in structures.get_all_structure_names():
        line = '<tr> <td>{0}</td> <td>{1}</td> <td>{2}</td> </tr>\n' \
            .format(i[3].decode('utf-8', errors='ignore'),
                    i[4].decode('utf-8', errors='ignore'),
                    i[2].decode('utf-8', errors='ignore'),
                    # i[0]
                    )
        # i[0] -> id
        table_string += line
        table_string = table_string
    p = pathlib.Path("./cgi/strflog_Template.htm")
    t = Template(p.read_bytes().decode('utf-8', 'ignore'))
    return str(t.safe_substitute({"logtablecolumns": table_string})).encode('utf-8', 'ignore')


if __name__ == "__main__":
    try:
        import wsgiref.simple_server

        server = wsgiref.simple_server.make_server('127.0.0.1', 8000, application)
        server.serve_forever()
        print("Webserver running...")
    except KeyboardInterrupt:
        print("Webserver stopped...")
