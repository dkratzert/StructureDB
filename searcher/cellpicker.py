"""
Created on 09.02.2015

 ----------------------------------------------------------------------------
* "THE BEER-WARE LICENSE" (Revision 42):
* <daniel.kratzert@uni-freiburg.de> wrote this file. As long as you retain this 
* notice you can do whatever you want with this stuff. If we meet some day, and 
* you think this stuff is worth it, you can buy me a beer in return.
* ----------------------------------------------------------------------------

@author: daniel
"""
import codecs
import os
import time
from pprint import pprint

import sys

from searcher import misc
from searcher.misc import get_error_from_value

essential_fields = ('_cell_length_a', '_cell_length_b', '_cell_length_c', '_cell_angle_alpha', '_cell_angle_beta',
                         '_cell_angle_gamma')
fields = ('_diffrn_ambient_temperature', '_diffrn_reflns_av_R_equivalents', '_diffrn_reflns_av_unetI/netI',
               '_diffrn_reflns_theta_max')


def get_cif_cell_raw(filename):
    """
    Get the cell of a cif file
    _cell_length_a                    29.227(4)
    _cell_length_b                    6.6568(8)
    _cell_length_c                    11.8204(14)
    _cell_angle_alpha                 90
    _cell_angle_beta                  107.055(5)
    _cell_angle_gamma                 90
    and the respective error values
    ...
    _shelx_hkl_file
    ;
    hkl data
    ;
    Attention: This implementation currently uses only the first cell of a cif file
    
    - list of desired fileds
    - cifdict = {}
    for line in f:
        for f in fields:
            test = line[:len(f)]
            if test == f:
                cifdict[f] = test 
    
    """
    name, a, b, c, alpha, beta, gamma, esda, esdb, esdc, esdalpha, esdbeta, esdgamma = \
        None, None, None, None, None, None, None, None, None, None, None, None, None
    filename = os.path.abspath(filename.name)
    cell = [None,  # 0 name
            None,  # 1 a
            None,  # 2 b
            None,  # 3 c
            None,  # 4 alpha
            None,  # 5 beta
            None,  # 6 gamma
            None,  # 7 esda
            None,  # 8 esdb
            None,  # 9 esdc
            None,  # 10 esdalpha
            None,  # 11 esdbeta
            None]  # 12 esdgamma
    with codecs.open(filename, "r", encoding='ascii', errors='ignore') as f:
        for line in f:
            if line.startswith('data_'):
                name = line.split('_')[1].strip('\n')
                cell[0] = name
            if line.startswith('_cell_length_a'):
                a = line.split()[1].split('(')[0]
                cell[1] = float(a)
                cell[7] = get_error_from_value(a)
            if line.startswith('_cell_length_b'):
                b = line.split()[1].split('(')[0]
                cell[2] = float(b)
                cell[8] = get_error_from_value(b)
            if line.startswith('_cell_length_c'):
                c = line.split()[1].split('(')[0]
                cell[3] = float(c)
                cell[9] = get_error_from_value(c)
            if line.startswith('_cell_angle_alpha'):
                alpha = line.split()[1].split('(')[0]
                cell[4] = float(alpha)
                cell[10] = get_error_from_value(alpha)
            if line.startswith('_cell_angle_beta'):
                beta = line.split()[1].split('(')[0]
                cell[5] = float(beta)
                cell[11] = get_error_from_value(beta)
            if line.startswith('_cell_angle_gamma'):
                gamma = line.split()[1].split('(')[0]
                cell[6] = float(gamma)
                cell[12] = get_error_from_value(gamma)
            #if line.startswith("_shelx_hkl_file"):
            #    return cell
            if all(cell[n] for n, i in enumerate(cell)):
                # contains 1+n cells. n=0-inf
                return cell
    return []


class Cif():
    def __init__(self, file):
        """
        A cif file parsing object optimized for speed and simplicity.
        :param file: input filename object
        :type file: object
        
        - find loops
          - get loop header
          - start and end
        - fill for example atom loop dict with atoms
        """
        self.cif_data = {}
        self.essential_fields = essential_fields
        self.fields = fields
        self.all_fields = self.essential_fields + self.fields
        self.ok = self.parsefile(file)

    def parsefile(self, file):
        """
        This method parses the cif file. Currently, only single items and atoms are supported.
        
        :param file: Cif file name
        :type file: str
        :return: cif file content
        :rtype: dict
        """
        data = False
        loophead_list = []
        loop = False
        atoms = {}
        atkey = ''
        with open(file, mode='r', encoding='ascii', errors="ignore") as f:
            for num, line in enumerate(f):
                line = line.lstrip().strip('\r\n ')
                if loop:
                    # Do not parse atom type symbols:
                    if line == "_atom_type_symbol":
                        loop = False
                        continue
                    if line == '_atom_site_label':
                        atkey = '_atom_site_label'
                    if line == '_atom_site_aniso_label':
                        atkey = '_atom_site_aniso_label'
                    if not line:
                        loop = False
                        loophead_list.clear()
                        atkey = ''
                        continue
                    # collecting keywords from head:
                    if line[:1] == "_":
                        loophead_list.append(line)
                    # We are in a loop and the header ended, so we collect data:
                    else:
                        loopitem = {}
                        loop_data_line = line.split()
                        atom = loop_data_line[0]
                        # quick hack, have to unwrap wrapped loop data:
                        if len(loop_data_line) != len(loophead_list):
                            # parse lines until next _ at line start
                            continue
                        for n, item in enumerate(loop_data_line):
                            loopitem[loophead_list[n]] = item
                        if atkey and loopitem[atkey] in atoms:
                            # atom is already there not there, upating values
                            atoms[loopitem[atkey]].update(loopitem)
                        elif atkey:
                            # atom is not there, creating key
                            atoms[loopitem[atkey]] = loopitem
                # First find the start of the cif (the data tag)
                if line[:5] == 'data_':
                    if not data:
                        name = line.split('_')[1].strip('\n\r')
                        self.cif_data['data'] = name
                        data = True
                        continue
                    else:
                        break  # TODO: support multi cif
                # Find the loop positions:
                if line[:5] == "loop_":
                    loop = True
                    continue
                if line.startswith('_') and not loop:
                    lsplit = line.split()
                    if len(lsplit) > 1:
                        self.cif_data[lsplit[0]] = lsplit[1]
                if line[:15] == "_shelx_hkl_file":
                    break
                if line[:7] == "_refln_":
                    break
        self.cif_data['_atom'] = atoms
        #pprint(self.cif_data)
        if not data:
            return False
        if not atoms:
            self.cif_data.clear()
            return False
        else:
            return True

    def __iter__(self):
        """
        an iterable for the Cif object
        :return: cif entries
        """
        if self.ok:
            yield self.cif_data
        else:
            yield {}

    @property
    def _cell_length_a(self):
        if "_cell_length_a" in self.cif_data:
            return self.cif_data["_cell_length_a"]
        else:
            return ''

    @property
    def _cell_length_b(self):
        if "_cell_length_b" in self.cif_data:
            return self.cif_data["_cell_length_b"]
        else:
            return ''

    @property
    def _cell_length_c(self):
        if "_cell_length_c" in self.cif_data:
            return self.cif_data["_cell_length_c"]
        else:
            return ''

    @property
    def _cell_angle_alpha(self):
        if "_cell_length_alpha" in self.cif_data:
            return self.cif_data["_cell_length_alpha"]
        else:
            return ''

    @property
    def _cell_angle_beta(self):
        if "_cell_length_beta" in self.cif_data:
            return self.cif_data["_cell_length_beta"]
        else:
            return ''

    @property
    def _cell_angle_gamma(self):
        if "_cell_length_gamma" in self.cif_data:
            return self.cif_data["_cell_length_gamma"]
        else:
            return ''


def get_res_cell(filename):
    """
    Returns the unit cell parameters from the list file as list:
    ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
    """
    file_list = misc.open_file_read(filename)
    cell = False
    for line in file_list:
        if line[:4] == 'CELL':
            cell = line.split()[2:8]
            try:
                cell = [float(i) for i in cell]
                if not len(cell) == 6:
                    raise ValueError
            except ValueError:
                print('Bad cell parameters in {0}.'.format(filename))
                return False
            if line[:4] == 'UNIT':
                break
    if not cell:
        #print('Unable to find unit cell parameters in the file.')
        return False
    return cell


if __name__ == '__main__':
    time1 = time.clock()
    #c = Cif("test-data/p21c.cif")
    #import CifFile
    #c = CifFile.ReadCif("test-data/p21c.cif")
    cifffile = "/Users/daniel/.olex2/data/3e30b45376c2d4175951f811f7137870/customisation/cif_templates/ALS_BL1131_post_07_2014.cif"
    cif2 = "/Users/daniel/Documents/Strukturen/DK_ML7-66/DK_ML7-66-final-old.cif"
    cif3 = "/Users/daniel/Downloads/.olex/originals/10000007.cif"
    c = Cif(cif3)
    time2 = time.clock()
    diff = round(time2-time1, 4)
    print(diff, 's')
    #sys.exit()
    for i in c:
        pprint(i)
