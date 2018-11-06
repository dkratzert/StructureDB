#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on 09.02.2015

 ----------------------------------------------------------------------------
* "THE BEER-WARE LICENSE" (Revision 42):
* <daniel.kratzert@uni-freiburg.de> wrote this file. As long as you retain this 
* notice you can do whatever you want with this stuff. If we meet some day, and 
* you think this stuff is worth it, you can buy me a beer in return.
* ----------------------------------------------------------------------------

@author: Daniel Kratzert
"""
import fnmatch
import os
import pathlib
import re
import sys
import tarfile
import time
import zipfile

from searcher import atoms, database_handler
from lattice.lattice import vol_unitcell
from searcher.fileparser import Cif
from shelxfile.shelx import ShelXFile
DEBUG = False

excluded_names = ['ROOT',
                  '.OLEX',
                  'olex',
                  'TMP',
                  'TEMP',
                  'Papierkorb',
                  'Recycle.Bin',
                  'dsrsaves',
                  'BrukerShelXlesaves',
                  'shelXlesaves'
                  ]


class MyZipBase(object):
    def __init__(self, filepath: str) -> None:
        self.filepath = filepath
        self.cifname = ''
        self.cifpath = ''


class MyZipReader(MyZipBase):
    def __init__(self, filepath):
        """
        extracts .cif files from zip files
        """
        super().__init__(filepath)

    def __iter__(self) -> list:
        """
        returns an iterator of cif files in the zipfile as list.
        """
        try:
            zfile = zipfile.ZipFile(self.filepath)
            for name in zfile.namelist():
                (self.cifpath, self.cifname) = os.path.split(name)
                if self.cifname.endswith('.cif'):
                    if not self.cifname.startswith('__') and zfile.NameToInfo[name].file_size < 150000000:
                        yield zfile.read(name).decode('utf-8', 'ignore').splitlines(keepends=True)
        except Exception as e:
            #print("Error: '{}' in file {}".format(e, self.filepath.encode(encoding='utf-8', errors='ignore')))
            #print(e, self.filepath)  # filepath is not utf-8 save
            yield []


class MyTarReader(MyZipBase):
    def __init__(self, filepath):
        """
        extracts .cif files from tar.gz files
        """
        super().__init__(filepath)

    def __iter__(self) -> list:
        """
        returns an iterator of cif files in the zipfile as list.
        """
        try:
            tfile = tarfile.open(self.filepath, mode='r')
            for name in tfile.getnames():
                self.cifpath, self.cifname = os.path.split(name)
                if self.cifname.endswith('.cif'):
                    yield tfile.extractfile(name).read().decode('utf-8', 'ignore').splitlines(keepends=True)
        except Exception as e:
            #print("Error: '{}' in file {}".format(e, self.filepath.encode(encoding='utf-8', errors='ignore')))
            #print(e, self.filepath)  # filepath is not utf-8 save
            yield []


def create_file_list(searchpath='None', ending='cif'):
    """
    walks through the file system and collects cells from res/cif files.
    Pathlib is nice, but does not allow me to do rglob for more than one file type.
    """
    if not os.path.isdir(searchpath):
        print('search path {0} not found! Or no directory!'.format(searchpath))
        sys.exit()
    print('collecting files... (may take some minutes)')
    p = pathlib.Path(searchpath)
    paths = p.rglob("*.{}".format(ending))
    return paths


def filewalker_walk(startdir: str, patterns: list):
    """
    walks through the filesystem starting from startdir and searches
    for files with ending endings.

    Since os.walk() uses scandir, it is as fast as pathlib.
    """
    filelist = []
    print('collecting files below ' + startdir)
    for root, _, files in os.walk(startdir):
        for filen in files:
            omit = False
            if any(fnmatch.fnmatch(filen, pattern) for pattern in patterns):
                fullpath = os.path.abspath(os.path.join(root, filen))
                if os.stat(fullpath).st_size == 0:
                    continue
                for ex in excluded_names:
                    if re.search(ex, fullpath, re.I):
                        omit = True
                if omit:
                    continue
                if filen == 'xd_geo.cif':  # Exclude xdgeom cif files
                    continue
                if filen == 'xd_four.cif':  # Exclude xdfourier cif files
                    continue
                # This is much faster than yield():
                filelist.append([root, filen])
            else:
                continue
    return filelist


def put_files_in_db(self=None, searchpath: str = './', excludes: list = None, lastid: int = 1,
                    structures=None, fillcif=True, fillres=True) -> int:
    """
    Imports files from a certain directory
    :param fillres: Should it index res files or not.
    :param fillcif: Should it index cif files or not.
    """
    if excludes:
        excluded_names.extend(excludes)
    if not searchpath:
        return 0
    if self:
        structures = self.structures
    if lastid <= 1:
        lastid = 1
    prognum = 0
    num = 1
    zipcifs = 0
    rescount = 0
    cifcount = 0
    time1 = time.clock()
    patterns = ['*.cif', '*.zip', '*.tar.gz', '*.tar.bz2', '*.tgz', '*.res']
    filelist = filewalker_walk(str(searchpath), patterns)
    options = {}
    filecount = 1
    for filenum, (filepth, name) in enumerate(filelist, start=1):
        filecount = filenum
        fullpath = os.path.join(filepth, name)
        options['modification_time'] = time.strftime('%Y-%m-%d', time.gmtime(os.path.getmtime(fullpath)))
        options['file_size'] = int(os.stat(str(fullpath)).st_size)
        cif = Cif(options=options)
        if self:
            if prognum == 20:
                prognum = 0
            self.progressbar(prognum, 0, 20)
        # This is really ugly copy&pase code. TODO: refractor this:
        if name.endswith('.cif') and fillcif:
            with open(fullpath, mode='r', encoding='ascii', errors="ignore") as f:
                try:
                    cifok = cif.parsefile(f.readlines())
                    if not cifok:
                        if DEBUG:
                            print("Could not parse: {}.".format(fullpath.encode('ascii', 'ignore')))
                        continue
                except IndexError:
                    continue
                if cif:  # means cif object has data inside (cif could be parsed)
                    tst = fill_db_tables(cif, filename=name, path=filepth, structure_id=lastid,
                                         structures=structures)
                    if not tst:
                        continue
                    if self:
                        self.add_table_row(name, filepth, cif.cif_data['data'], str(lastid))
                    cifcount += 1
                    lastid += 1
                    num += 1
                    if lastid % 1000 == 0:
                        print('{} files ...'.format(num))
                        structures.database.commit_db()
                    prognum += 1
            continue
        if (name.endswith('.zip') or name.endswith('.tar.gz') or name.endswith('.tar.bz2')
                or name.endswith('.tgz')) and fillcif:
            if fullpath.endswith('.zip'):
                # MyZipReader defines .cif ending:
                z = MyZipReader(fullpath)
            else:
                z = MyTarReader(fullpath)
            for zippedfile in z:              # the list of cif files in the zip file
                omit = False
                for ex in excluded_names:          # remove excludes
                    if re.search(ex, z.cifpath, re.I):
                        omit = True
                if omit:
                    continue
                try:
                    cifok = cif.parsefile(zippedfile)
                    if not cifok:
                        if DEBUG:
                            print("Could not parse: {}.".format(fullpath.encode('ascii', 'ignore')))
                        continue
                except IndexError:
                    continue
                if cif:
                    tst = fill_db_tables(cif, filename=z.cifname, path=fullpath,
                                         structure_id=str(lastid), structures=structures)
                    if not tst:
                        if DEBUG:
                            print('cif file not added:', fullpath) 
                        continue
                    if self:
                        self.add_table_row(name=z.cifname, path=fullpath,
                                           data=cif.cif_data['data'], structure_id=str(lastid))
                    zipcifs += 1
                    cifcount += 1
                    lastid += 1
                    num += 1
                    if lastid % 1000 == 0:
                        print('{} files ...'.format(num))
                        structures.database.commit_db()
                    prognum += 1
            continue
        if name.endswith('.res') and fillres:
            tst = None
            try:
                res = ShelXFile(fullpath)
            except Exception as e:
                if DEBUG:
                    print(e)
                    print("Could not parse: {}.".format(fullpath.encode('ascii', 'ignore')))
                continue
            if res:
                tst = fill_db_with_res_data(res, filename=name, path=filepth, structure_id=lastid,
                                            structures=structures, options=options)
            if not tst:
                if DEBUG:
                    print('res file not added:', fullpath)
                continue
            if self:
                self.add_table_row(name=name, path=fullpath, data=name, structure_id=str(lastid))
            lastid += 1
            num += 1
            rescount += 1
            if lastid % 1000 == 0:
                print('{} files ...'.format(num))
                structures.database.commit_db()
            prognum += 1
            continue
        if self:
            if not self.decide_import:
                # This means, import was aborted.
                self.abort_import_button.hide()
                self.decide_import = True
                break
    structures.database.commit_db()
    time2 = time.clock()
    diff = time2 - time1
    m, s = divmod(diff, 60)
    h, m = divmod(m, 60)
    tmessage = 'Added {0} files ({5} cif, {6} res) files ({4} in compressed files) to database in: ' \
               '{1:>2d} h, {2:>2d} m, {3:>3.2f} s'
    print('      {} files considered.'.format(filecount))
    print(tmessage.format(num - 1, int(h), int(m), s, zipcifs, cifcount, rescount))
    if self:
        self.ui.statusbar.showMessage(tmessage.format(num - 1, int(h), int(m), s, zipcifs, cifcount, rescount))
    return lastid-1


def fill_db_tables(cif: Cif, filename: str, path: str, structure_id: str,
                   structures: database_handler.StructureTable):
    """
    Fill all info from cif file into the database tables
    _atom_site_label
    _atom_site_type_symbol
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    _atom_site_U_iso_or_equiv
    _atom_site_adp_type
    _atom_site_occupancy
    _atom_site_site_symmetry_order
    _atom_site_calc_flag
    _atom_site_refinement_flags_posn
    _atom_site_refinement_flags_adp
    _atom_site_refinement_flags_occupancy
    _atom_site_disorder_assembly
    _atom_site_disorder_group
    """
    a = cif._cell_length_a
    b = cif._cell_length_b
    c = cif._cell_length_c
    alpha = cif._cell_angle_alpha
    beta = cif._cell_angle_beta
    gamma = cif._cell_angle_gamma
    volume = cif._cell_volume
    if not all((a, b, c, alpha, beta, gamma)):
        return False
    if not volume or volume == "?":
        # TODO: bring get_error_from_value() to here:
        try:
            if isinstance(a, str):
                a = float(a.split('(')[0])
            if isinstance(b, str):
                b = float(b.split('(')[0])
            if isinstance(c, str):
                c = float(c.split('(')[0])
            if isinstance(alpha, str):
                alpha = float(alpha.split('(')[0])
            if isinstance(beta, str):
                beta = float(beta.split('(')[0])
            if isinstance(gamma, str):
                gamma = float(gamma.split('(')[0])
            volume = str(vol_unitcell(a, b, c, alpha, beta, gamma))
        except ValueError:
            volume = ''
    measurement_id = structures.fill_measuremnts_table(filename, structure_id)
    structures.fill_structures_table(path, filename, structure_id, measurement_id, cif.cif_data['data'])
    structures.fill_cell_table(structure_id, a, b, c, alpha, beta, gamma, volume)
    sum_from_dict = {}
    for x in cif._atom:
        try:
            try:
                disord = int(cif._atom[x]['_atom_site_disorder_group'])
            except (KeyError, ValueError):
                disord = 0
            try:
                occu = float(cif._atom[x]['_atom_site_occupancy'].split('(')[0])
            except (KeyError, ValueError):
                occu = 1.0
            try:
                atom_type_symbol = atoms.get_atomlabel(cif._atom[x]['_atom_site_type_symbol'])
            except KeyError:
                atom_type_symbol  = atoms.get_atomlabel(x)
            elem = atom_type_symbol.capitalize()
            structures.fill_atoms_table(structure_id, x,
                                         atom_type_symbol,
                                         cif._atom[x]['_atom_site_fract_x'].split('(')[0],
                                         cif._atom[x]['_atom_site_fract_y'].split('(')[0],
                                         cif._atom[x]['_atom_site_fract_z'].split('(')[0],
                                         occu,
                                         disord
                                        )
            if elem in sum_from_dict:
                sum_from_dict[elem] += occu
            else:
                sum_from_dict[elem] = occu
        except KeyError as e:
            #print(x, filename, e)
            pass
    cif.cif_data['calculated_formula_sum'] = sum_from_dict
    structures.fill_residuals_table(structure_id, cif)
    return True


def fill_db_with_res_data(res: ShelXFile, filename: str, path: str, structure_id: str,
                          structures: database_handler.StructureTable, options: dict):
    if not res.cell:
        return False
    if not all([res.cell.a, res.cell.b, res.cell.al, res.cell.be, res.cell.ga]):
        return False
    if not res.cell.volume:
        return False
    measurement_id = structures.fill_measuremnts_table(filename, structure_id)
    structures.fill_structures_table(path, filename, structure_id, measurement_id, res.titl)
    structures.fill_cell_table(structure_id, res.cell.a, res.cell.b, res.cell.c, res.cell.al,
                               res.cell.be, res.cell.ga, res.cell.volume)
    for at in res.atoms:
        if at.qpeak:
            continue
        if at.element.lower() == 'cnt':  # Do not add Shelxle centroids
            continue
        structures.fill_atoms_table(structure_id, 
                                    at.name,
                                    at.element.capitalize(),
                                    at.x,
                                    at.y,
                                    at.z,
                                    at.sof,
                                    at.part.n)
    cif = Cif(options=options)
    cif.cif_data["_cell_formula_units_Z"] = res.Z
    cif.cif_data["_space_group_symop_operation_xyz"] = "\n".join([repr(x) for x in res.symmcards])
    try:
        cif.cif_data["calculated_formula_sum"] = res.sum_formula_ex_dict()
    except ZeroDivisionError:
        pass
    try:
        cif.cif_data["_chemical_formula_sum"] = res.sum_formula_exact
    except ZeroDivisionError:
        pass
    cif.cif_data["_diffrn_radiation_wavelength"] = res.wavelen
    if res.R1:
        cif.cif_data["_refine_ls_R_factor_gt"] = res.R1
    if res.wR2:
        cif.cif_data["_refine_ls_wR_factor_ref"] = res.wR2
    if res.parameters:
        cif.cif_data['_refine_ls_number_parameters'] = res.parameters
    if res.data:
        cif.cif_data['_refine_ls_number_reflns'] = res.data
    if res.num_restraints:
        cif.cif_data['_refine_ls_number_restraints'] = res.num_restraints
    if res.temp_in_Kelvin:
        cif.cif_data['_diffrn_ambient_temperature'] = round(res.temp_in_Kelvin, 5)
    if res.dhole:
        cif.cif_data['_refine_diff_density_min'] = res.dhole
    if res.hpeak:
        cif.cif_data['_refine_diff_density_max'] = res.hpeak
    if res.latt:
        cif.cif_data['_space_group_centring_type'] = res.latt.N_str
    if res.space_group:
        cif.cif_data["_space_group_name_H-M_alt"] = res.space_group
    if res.goof:
        cif.cif_data["_refine_ls_goodness_of_fit_ref"] = res.goof
    if res.rgoof:
        cif.cif_data["_refine_ls_restrained_S_all"] = res.rgoof
    try:
        cif.cif_data["_shelx_res_file"] = str(res)
    except IndexError:
        pass
    structures.fill_residuals_table(structure_id, cif)
    return True


if __name__ == '__main__':
    z = MyTarReader('./test-data/106c.tar.bz2')

    for i in z:
        print(i)

    #filewalker_walk('./')
    #z = zipopener('../test-data/Archiv.zip')
    #print(z)

    #fp = create_file_list('../test-data/', 'zip')
    #for i in fp:
    #    print(i)