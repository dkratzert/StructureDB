# -*- coding: utf-8 -*-
# !C:\tools\Python-3.6.2_64\pythonw.exe
# !/usr/local/bin/python3.6

###########################################################
###  Configure the web server here:   #####################
from contextlib import contextmanager
from os.path import isfile

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from searcher.database_handler import get_residuals, get_cell_by_id, find_by_it_number, find_by_strings, \
    get_atoms_table, find_by_elements, find_by_date, find_cell_by_volume, get_cells_as_list

host = "10.6.13.3"
port = "80"
dbfilename = "../structurefinder.sqlite"

###########################################################

import socket
names = ['PC9', 'DDT-2.local']
# run on local ip on my PC:
if socket.gethostname() in names:
    host = '127.0.0.1'
    port = "8080"
site_ip = host + ':' + port

import math
import os
import pathlib
import sys


pyver = sys.version_info
if pyver[0] == 3 and pyver[1] < 4:
    # Python 2 creates a syntax error anyway.
    print("You need Python 3.4 and up in oder to run this proram!")
    sys.exit()

try:  # Adding local path to PATH
    sys.path.insert(0, os.path.abspath('./'))
except(KeyError, ValueError):
    print('Unable to set PATH properly. strf_web.py might not work.')

from cgi_ui.bottle import Bottle, static_file, template, redirect, request, response
from displaymol.mol_file_writer import MolFile
from displaymol.sdm import SDM
from lattice import lattice
from pymatgen.core import mat_lattice
from searcher.misc import is_valid_cell, flatten, is_a_nonzero_file, get_list_of_elements

"""
TODO:
- Make login infrastructure.
- Add option: should contain *only* these elements
- Maybe http://www.daterangepicker.com
"""

if isfile(dbfilename):
    engine = create_engine('sqlite:///' + dbfilename)
    Session = sessionmaker(bind=engine)


@contextmanager
def session_scope():
    """Provide a transactional scope around a series of operations."""
    session = Session()
    # session.autoflush = False
    try:
        yield session
        session.commit()
    except:
        session.rollback()
        raise
    finally:
        session.close()


app = application = Bottle()
# bottle.debug(True)  # Do not enable debug in production systems!


@app.route('/all')
def structures_list_data():
    """
    The content of the structures list.
    """
    with session_scope() as session:
        return get_structures_json(session, show_all=True)


@app.route('/', method=['POST', 'GET'])
def main():
    """
    The main web site with html template and space group listing.
    """
    response.content_type = 'text/html; charset=UTF-8'
    p = pathlib.Path('./cgi_ui/views/spgr.html').open()
    space_groups = p.read().encode(encoding='UTF-8', errors='ignore')
    p.close()
    output = template('./cgi_ui/views/strf_web_template', {"my_ip": site_ip, "space_groups": space_groups})
    return output


@app.route("/cellsrch")
def cellsrch():
    cell_search = request.GET.cell_search
    more_results = (request.GET.more == "true")
    sublattice = (request.GET.supercell == "true")
    cell = is_valid_cell(cell_search)
    print("Cell search:", cell)
    if cell:
        with session_scope() as session:
            ids = find_cell(session, cell, more_results=more_results, sublattice=sublattice)
            print("--> Got {} structures from cell search.".format(len(ids)))
            return get_structures_json(session, ids, show_all=False)


@app.route("/txtsrch")
def txtsrch():
    text_search = request.GET.text_search
    print("Text search:", text_search)
    with session_scope() as session:
        ids = search_text(session, text_search)
        return get_structures_json(session, ids, show_all=False)


@app.route("/adv_srch")
def adv():
    elincl = request.GET.elements_in
    elexcl = request.GET.elements_out
    date1 = request.GET.date1
    date2 = request.GET.date2
    cell_search = request.GET.cell_search
    txt_in = request.GET.text_in
    txt_out = request.GET.text_out
    more_results = (request.GET.more == "true")
    sublattice = (request.GET.supercell == "true")
    it_num = request.GET.it_num
    print("Advanced search:", elincl, elexcl, date1, date2, cell_search, txt_in, txt_out, more_results, sublattice,
          it_num)
    with session_scope() as session:
        ids = advanced_search(session, cellstr=cell_search, elincl=elincl, elexcl=elexcl, txt_in=txt_in, txt_out=txt_out,
                              sublattice=sublattice, more_results=more_results, date1=date1, date2=date2, it_num=it_num)
        print("--> Got {} structures from Advanced search.".format(len(ids)))
        return get_structures_json(session, ids)


@app.route('/molecule', method='POST')
def jsmol_request():
    """
    A request for atom data from jsmol.
    """
    str_id = request.POST.id
    print("Molecule id:", str_id)
    with session_scope() as session:
        if str_id:
            cell = get_cell_by_id(session, str_id)
            if request.POST.grow == 'true':
                symmcards = [x.split(',') for x in get_row_as_dict(str_id)
                ['_space_group_symop_operation_xyz'].replace("'", "").replace(" ", "").split("\n")]
                atoms = get_atoms_table(str_id, cell[:6], cartesian=False)
                if atoms:
                    sdm = SDM(atoms, symmcards, cell)
                    needsymm = sdm.calc_sdm()
                    atoms = sdm.packer(sdm, needsymm)
            else:
                atoms = get_atoms_table(session, str_id, cell[:6], cartesian=True)
            try:
                m = MolFile(atoms)
                return m.make_mol()
            except(KeyError, TypeError) as e:
                print('Exception in jsmol_request: {}'.format(e))
                return ''


@app.route('/', method='POST')
def post_request():
    """
    Handle POST requests.
    """
    cif_dic = {}
    str_id = request.POST.id
    resid1 = request.POST.residuals1 == 'true'
    resid2 = request.POST.residuals2 == 'true'
    all_cif = (request.POST.all == 'true')
    unitcell = request.POST.unitcell
    print("Structure id:", str_id)
    with session_scope() as session:
        if str_id:
            cif_dic = get_residuals(session, str_id)
        if str_id and unitcell and not (resid1 or resid2 or all_cif):
            try:
                return get_cell_by_id(session, str_id)
            except ValueError as e:
                print("Exception raised:")
                print(e)
                return ''
        if str_id and resid1:
            return get_residuals_table1(cif_dic)
        if str_id and resid2:
            return get_residuals_table2(cif_dic)
        if str_id and all_cif:
            return get_all_cif_val_table(session, str_id)


@app.route('/static/<filepath:path>')
def server_static(filepath):
    """
    Static files such as images or CSS files are not served automatically.
    The static_file() function is a helper to serve files in a safe and convenient way (see Static Files).
    This example is limited to files directly within the /path/to/your/static/files directory because the
    <filename> wildcard won’t match a path with a slash in it. To serve files in subdirectories, change
    the wildcard to use the path filter:
    """
    return static_file(filepath, root='./cgi_ui/static/')


@app.route('/version')
def version():
    from misc.version import VERSION
    return 'version ' + str(VERSION)


@app.route('/cgi-bin/strf_web.cgi')
def redirect_old_path():
    redirect('/')


@app.route('/favicon.ico')
def redirect_to_favicon():
    redirect('/static/favicon.ico')


@app.error(404)
def error404(error):
    """
    Redefine 404 message.
    """
    return '''<div style="text-align: center;">
                <b>Nothing here, sorry.</b><br>
                <p>
                <a href="http://{}{}/">Back to main page</a>
                </p>
              </div>
            '''.format(host, ':' + port)


def is_ajax():
    if request.headers.get('X-Requested-With') == 'XMLHttpRequest':
        return True
    else:
        return False


def get_structures_json(session: Session, ids: (list, tuple) = None, show_all: bool = False) -> dict:
    """
    Returns the next package of table rows for continuos scrolling.
    """
    failure = {
        "status" : "error",
        "message": "Nothing found."
    }
    if not ids and not show_all:
        # return json.dumps(failure)
        return {}
    dic = get_all_structures_as_dict(ids, all_ids=show_all)
    number = len(dic)
    print("--> Got {} structures from actual search.".format(number))
    if number == 0:
        # return json.dumps(failure)
        return {}
    return {"total": number, "records": dic, "status": "success"}


def get_cell_parameters(session: Session, strid: str) -> str:
    """
    Resturns unit cell parameters as html formated string.
    """
    c = get_cell_by_id(session, strid)
    cstr = """<b>Unit Cell:</b>&nbsp;&nbsp; 
                      <i>a</i> = {0:>8.3f}&nbsp;&angst;,&nbsp;
                      <i>b</i> = {1:>8.3f}&nbsp;&angst;,&nbsp;
                      <i>c</i> = {2:>8.3f}&nbsp;&angst;,&nbsp; 
                      <i>&alpha;</i> = {3:>8.3f}&deg;,&nbsp;
                      <i>&beta;</i> = {4:>8.3f}&deg;,&nbsp;
                      <i>&gamma;</i> = {5:>8.3f}&deg;,&nbsp;
                      <i>V</i> = {6}&nbsp;&angst;<sup>3</sup>&nbsp;&nbsp;&nbsp;&nbsp; 
            <div style="font-size:0pt" id='hidden-cell'>{0}  {1}  {2}  {3}  {4}  {5}</div>
            """.format(c[0], c[1], c[2], c[3], c[4], c[5], round(c[6], 2))
    return cstr


def get_residuals_table1(cif_dic: dict) -> str:
    """
    Returns a table with the most important residuals of a structure.
    """
    try:
        rsigma = " / {}".format(cif_dic['_diffrn_reflns_av_unetI_netI'])
    except (TypeError, ValueError):
        rsigma = " "
    if not cif_dic:
        return ""
    if cif_dic['_refine_diff_density_max']:
        peakhole = "{} / {}".format(cif_dic['_refine_diff_density_max'], cif_dic['_refine_diff_density_min'])
    else:
        peakhole = " "
    table1 = """
    <table class="table table-bordered" id='resitable1'>
        <tbody>
        <tr><td style='width: 40%'><b>Space Group</b></td>                 <td>{0}</td></tr>
        <tr><td><b>Z</b></td>                           <td>{1}</td></tr>
        <tr><td><b>Sum Formula</b></td>                 <td>{2}</td></tr>
        <tr><td><b>Temperature [K]</b></td>             <td>{3}</td></tr>
        <tr><td><b><i>wR</i><sub>2</sub></b></td>       <td>{4}</td></tr>
        <tr><td><b><i>R<i/><sub>1</sub></b></td>        <td>{5}</td></tr>
        <tr><td><b>Goof</b></td>                        <td>{6}</td></tr>
        <tr><td><b>Max Shift/esd</b></td>               <td>{7}</td></tr>
        <tr><td><b>Peak / Hole [e&angst;<sup>&minus;3</sup>]</b></td>             <td>{8}</td></tr>
        <tr><td><b><i>R</i><sub>int</sub> / <i>R</i><sub>&sigma;</sub></b></b></td>    <td>{9}{10} </td></tr>
        <tr><td><b>Wavelength [&angst;]</b></td>                      <td>{11}</td></tr>
        </tbody>
    </table>
    """.format(cif_dic['_space_group_name_H_M_alt'],
               cif_dic['_cell_formula_units_Z'],
               cif_dic['_chemical_formula_sum'],
               cif_dic['_diffrn_ambient_temperature'],
               cif_dic['_refine_ls_wR_factor_ref'],
               cif_dic['_refine_ls_R_factor_gt'],
               cif_dic['_refine_ls_goodness_of_fit_ref'],
               cif_dic['_refine_ls_shift_su_max'],
               peakhole,
               cif_dic['_diffrn_reflns_av_R_equivalents'],
               rsigma,
               cif_dic['_diffrn_radiation_wavelength']
               )
    return table1


def get_residuals_table2(cif_dic: dict) -> str:
    """
    Returns a table with the most important residuals of a structure.
    """
    # cell = structures.get_cell_by_id(structure_id)
    if not cif_dic:
        return ""
    wavelen = cif_dic['_diffrn_radiation_wavelength']
    thetamax = cif_dic['_diffrn_reflns_theta_max']
    # d = lambda/2sin(theta):
    try:
        d = wavelen / (2 * math.sin(math.radians(thetamax)))
    except(ZeroDivisionError, TypeError):
        d = 0.0
    try:
        compl = cif_dic['_diffrn_measured_fraction_theta_max'] * 100
        if not compl:
            compl = 0.0
        if isinstance(compl, str):
            compl = 0.0
    except TypeError:
        compl = 0.0
    try:
        data_to_param = cif_dic['_refine_ls_number_reflns'] / cif_dic['_refine_ls_number_parameters']
    except TypeError:
        data_to_param = 0
    table2 = """
    <table class="table table-bordered" id='resitable2'>
        <tbody>
        <tr><td style='width: 40%'><b>Measured Refl.</b></td>       <td>{0}</td></tr>
        <tr><td><b>Independent Refl.</b></td>                       <td>{9}</td></tr>
        <tr><td><b>Data with [<i>I</i>>2&sigma;(<i>I</i>)] </b></td>    <td>{10}</td></tr>
        <tr><td><b>Parameters</b></td>                              <td>{1}</td></tr>
        <tr><td><b>data/param</b></td>                              <td>{2:<5.1f}</td></tr>
        <tr><td><b>Restraints</b></td>                              <td>{3}</td></tr>
        <tr><td><b>&theta;<sub>max</sub> [&deg;]</b></td>                    <td>{4}</td></tr>
        <tr><td><b>&theta;<sub>full</sub> [&deg;]</b></td>                   <td>{5}</td></tr>
        <tr><td><b>d [&angst;]</b></td>                             <td>{6:5.3f}</td></tr>
        <tr><td><b>completeness [%]</b></td>                            <td>{7:<5.1f}</td></tr>
        <tr><td><b>CCDC Number</b></td>                             <td>{8}</td></tr>
        </tbody>
    </table>
    """.format(cif_dic['_diffrn_reflns_number'],
               cif_dic['_refine_ls_number_parameters'],
               data_to_param,
               cif_dic['_refine_ls_number_restraints'],
               thetamax,
               cif_dic['_diffrn_reflns_theta_full'],
               d,
               compl,
               cif_dic['_database_code_depnum_ccdc_archive'],
               cif_dic['_refine_ls_number_reflns'],
               cif_dic['_reflns_number_gt']
               )
    return table2


def get_all_cif_val_table(session, structure_id: int) -> str:
    """
    Returns a html table with the residuals values of a structure.
    """
    # starting table header (the div is for css):
    # style="white-space: pre": preserves white space
    table_string = """<h4>All CIF values</h4>
                        <div id="myresidualtable">
                        <table class="table table-striped table-bordered" style="white-space: pre">
                            <thead>
                                <tr>
                                    <th> Item </th>
                                    <th> Value </th>
                                </tr>
                            </thead>
                        <tbody>"""
    # get the residuals of the cif file as a dictionary:
    dic = get_residuals(session, structure_id)
    if not dic:
        return ""
    # filling table with data rows:
    for key, value in dic.items():
        if key == "Id":
            continue
        if isinstance(value, str):
            value = ''.join([x.replace("\n", "<br>").rstrip('\r\n') for x in value])
        if key == '_shelx_res_file':
            # Adding an ID to make font monospace in the view:
            table_string += '''<tr>
                                 <td class="residual-{}"> {} </a></td> 
                                 <td id=resfile > {} </a></td> 
                               </tr> \n'''.format(structure_id, key, value)
        else:
            table_string += '''<tr>
                                <td class="residual-{}"> {} </a></td> 
                                <td> {} </a></td> 
                           </tr> \n'''.format(structure_id, key, value)
    # closing table:
    table_string += """ </tbody>
                        </table>
                        </div>"""
    return table_string


def chunks(l: list, n: int) -> list:
    """
    returns successive n-sized chunks from l.
    >>> l = [1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 'a', 'b', 'c', 'd', 'e', 'f']
    >>> chunks(l, 5)
    [[1, 2, 3, 4, 5], [6, 7, 8, 9, 0], ['a', 'b', 'c', 'd', 'e'], ['f']]
    >>> chunks(l, 1)
    [[1], [2], [3], [4], [5], [6], [7], [8], [9], [0], ['a'], ['b'], ['c'], ['d'], ['e'], ['f']]
    >>> chunks(l, 50)
    [[1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 'a', 'b', 'c', 'd', 'e', 'f']]
    """
    return [l[i:i + n] for i in range(0, len(l), n)]


def find_cell(session: Session, cell: list, sublattice=False, more_results=False) -> list:
    """
    Finds unit cells in db. Rsturns hits a a list of ids.
    """
    if more_results:
        # more results:
        vol_threshold = 0.09
        ltol = 0.2
        atol = 2
    else:
        # regular:
        vol_threshold = 0.03
        ltol = 0.06
        atol = 1
    volume = lattice.vol_unitcell(*cell)
    idlist = find_cell_by_volume(session, volume, vol_threshold)
    if sublattice:
        # sub- and superlattices:
        for v in [volume * x for x in [2.0, 3.0, 4.0, 6.0, 8.0, 10.0]]:
            # First a list of structures where the volume is similar:
            idlist.extend(find_cell_by_volume(session, v, vol_threshold))
        idlist = list(set(idlist))
        idlist.sort()
    idlist2 = []
    # Real lattice comparing in G6:
    if idlist:
        lattice1 = mat_lattice.Lattice.from_parameters_niggli_reduced(*cell)
        cells = []
        # SQLite can only handle 999 variables at once:
        for cids in chunks(idlist, 500):
            cells.extend(get_cells_as_list(cids))
        for num, cell_id in enumerate(idlist):
            try:
                lattice2 = mat_lattice.Lattice.from_parameters(
                        float(cells[num][2]),
                        float(cells[num][3]),
                        float(cells[num][4]),
                        float(cells[num][5]),
                        float(cells[num][6]),
                        float(cells[num][7]))
            except ValueError:
                continue
            mapping = lattice1.find_mapping(lattice2, ltol, atol, skip_rotation_matrix=True)
            if mapping:
                idlist2.append(cell_id)
    if idlist2:
        return idlist2
    else:
        return []


def search_text(session: Session, search_string: str) -> tuple:
    """
    searches db for given text
    """
    idlist = []
    if len(search_string) == 0:
        return ()
    if len(search_string) >= 2:
        if "*" not in search_string:
            search_string = "{}{}{}".format('*', search_string, '*')
    try:
        #  bad hack, should make this return ids like cell search
        idlist = tuple([x[0] for x in find_by_strings(session, search_string)])
    except AttributeError as e:
        print("Exception in search_text:")
        print(e)
    return idlist


def search_elements(session: Session, elements: str, anyresult: bool = False) -> list:
    """
    list(set(l).intersection(l2))
    """
    res = []
    try:
        formula = get_list_of_elements(elements)
    except KeyError:
        print('Element search error!')
        return []
    try:
        res = find_by_elements(session, formula, anyresult=anyresult)
    except AttributeError:
        print('Element search error!')
        pass
    return list(res)


def find_dates(session: Session, date1: str, date2: str) -> list:
    """
    Returns a list if id between date1 and date2
    """
    if not date1:
        date1 = '0000-01-01'
    if not date2:
        date2 = 'NOW'
    result = find_by_date(session, date1, date2)
    return result


def advanced_search(session: Session, cellstr: str, elincl, elexcl, txt_in, txt_out, sublattice, more_results,
                    date1: str = None, date2: str = None, it_num: str = None) -> list:
    """
    Combines all the search fields. Collects all includes, all excludes ad calculates
    the difference.
    """
    excl = []
    incl = []
    date_results = []
    results = []
    it_results = []
    cell = []
    if cellstr:
        cell = is_valid_cell(cellstr)
    if cell:
        cellres = find_cell(session, cell, sublattice=sublattice, more_results=more_results)
        incl.append(cellres)
    if elincl:
        incl.append(search_elements(session, elincl))
    if date1 != date2:
        date_results = find_dates(session, date1, date2)
    if it_num:
        try:
            it_results = find_by_it_number(session, int(it_num))
        except ValueError:
            pass
    if txt_in:
        if len(txt_in) >= 2 and "*" not in txt_in:
            txt_in = '*' + txt_in + '*'
        idlist = find_by_strings(session, txt_in)
        try:
            incl.append([i[0] for i in idlist])
        except(IndexError, KeyError):
            incl.append([idlist])  # only one result
    if elexcl:
        excl.append(search_elements(session, elexcl, anyresult=True))
    if txt_out:
        if len(txt_out) >= 2 and "*" not in txt_out:
            txt_out = '*' + txt_out + '*'
        idlist = find_by_strings(session, txt_out)
        try:
            excl.append([i[0] for i in idlist])
        except(IndexError, KeyError):
            excl.append([idlist])  # only one result
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
    if excl:
        # excl list should not be in the resukts at all
        try:
            return list(results - set(flatten(excl)))
        except TypeError:
            return []
    return list(results)


if __name__ == "__main__":
    print("Running on Python version {}".format(sys.version))
    if not is_a_nonzero_file(dbfilename):
        print("Unable to start!")
        print("The database file '{}' does not exist.".format(os.path.abspath(dbfilename)))
        sys.exit()
    print('### Running with database "{}" ###'.format(os.path.abspath(dbfilename)))
    # plain python wsgiref server:
    # app.run(host=host, port=port, reloader=True)
    # gunicorn server: Best used behind an nginx proxy server: http://docs.gunicorn.org/en/stable/deploy.html
    # you need "pip3 install gunicorn" to run this:
    # The current database interface allows only one worker (have to go to sqlalchemy!)
    app.run(host=host, port=port, reload=True, server='gunicorn', accesslog='-', errorlog='-', workers=1,
            access_log_format='%(h)s %(l)s %(u)s %(t)s "%(r)s" %(s)s %(b)s "%(f)s"')
