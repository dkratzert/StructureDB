# -*- coding: utf-8 -*-
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
import os
from pprint import pprint

from displaymol.sdm import SymmCards
from shelxfile.dsrmath import SymmetryElement


class Cif(object):
    def __init__(self, options=None):
        """
        A cif file parsing object optimized for speed and simplicity.
        It can not handle multi cif files.
        :param file: input filename object
        :type file: Path
        """
        if options is None:
            options = {'modification_time': "", 'file_size': ""}
        self.cif_data = {
            "data": '',
            "_cell_length_a": '',
            '_cell_length_b': '',
            '_cell_length_c': '',
            '_cell_angle_alpha': '',
            '_cell_angle_beta': '',
            '_cell_angle_gamma': '',
            "_cell_volume": '',
            "_cell_formula_units_Z": '',
            "_space_group_name_H-M_alt": '',
            "_space_group_name_Hall": '',
            "_space_group_centring_type": '',
            "_space_group_IT_number": '',
            "_space_group_crystal_system": '',
            "_space_group_symop_operation_xyz": '',
            "_audit_creation_method": '',
            "_chemical_formula_sum": '',
            "_chemical_formula_weight": '',
            "_exptl_crystal_description": '',
            "_exptl_crystal_colour": '',
            "_exptl_crystal_size_max": '',
            "_exptl_crystal_size_mid": '',
            "_exptl_crystal_size_min": '',
            "_exptl_absorpt_coefficient_mu": '',
            "_exptl_absorpt_correction_type": '',
            "_diffrn_ambient_temperature": '',
            "_diffrn_radiation_wavelength": '',
            "_diffrn_radiation_type": '',
            "_diffrn_source": '',
            "_diffrn_measurement_device_type": '',
            "_diffrn_reflns_number": '',
            "_diffrn_reflns_av_R_equivalents": '',
            "_diffrn_reflns_av_unetI/netI": '',
            "_diffrn_reflns_theta_min": '',
            "_diffrn_reflns_theta_max": '',
            "_diffrn_reflns_theta_full": '',
            "_diffrn_measured_fraction_theta_max": '',
            "_diffrn_measured_fraction_theta_full": '',
            "_reflns_number_total": '',
            "_reflns_number_gt": '',
            "_reflns_threshold_expression": '',
            "_reflns_Friedel_coverage": '',
            "_computing_structure_solution": '',
            "_computing_structure_refinement": '',
            "_refine_special_details": '',
            "_refine_ls_structure_factor_coef": '',
            "_refine_ls_weighting_details": '',
            "_refine_ls_number_reflns": '',
            "_refine_ls_number_parameters": '',
            "_refine_ls_number_restraints": '',
            "_refine_ls_R_factor_all": '',
            "_refine_ls_R_factor_gt": '',
            "_refine_ls_wR_factor_ref": '',
            "_refine_ls_wR_factor_gt": '',
            "_refine_ls_goodness_of_fit_ref": '',
            "_refine_ls_restrained_S_all": '',
            "_refine_ls_shift/su_max": '',
            "_refine_ls_shift/su_mean": '',
            "_refine_diff_density_max": '',
            "_refine_diff_density_min": '',
            "_diffrn_reflns_av_unetI_netI": '',
            "_database_code_depnum_ccdc_archive": '',
            "_shelx_res_file": '',
            "modification_time": options['modification_time'],
            "file_size": options['file_size']
            }

    def parsefile(self, txt):
        """
        This method parses the cif file. Currently, only single items and atoms are supported.
        :param file: Cif file name
        :type file: Path
        :return: cif file content
        :rtype: dict
        """
        data = False
        loop = False
        hkl = False
        loophead_list = []
        save_frame = False
        atoms = {}
        atkey = ''
        loop_body = False
        num = 0
        symm = False
        symmlist = []
        semi_colon_text_field = ''
        semi_colon_text_list = []
        cont = False  # continue to next line if True
        textlen = len(txt)
        for num, line in enumerate(txt):
            line = line.rstrip('\r\n ')
            if not line:
                loop = False
                loophead_list.clear()
                atkey = ''
                symm = False
                continue
            if line[0] == "_" and loop_body:
                loop = False
                loop_body = False
                symm = False
            if loop:
                line = line.lstrip()
                # leave out comments:
                if line[0] == '#':
                    continue
                # to collect the two parts of an atom loop (have to do it more general):
                if line == '_atom_site_label':
                    atkey = '_atom_site_label'
                if line == '_atom_site_aniso_label':
                    atkey = '_atom_site_aniso_label'
                if line == "_space_group_symop_operation_xyz" or line == '_symmetry_equiv_pos_as_xyz':
                    symm = '_space_group_symop_operation_xyz'
                    continue
                if line[:5] == "loop_":
                    loop = True
                    loop_body = False
                    loophead_list.clear()
                    atkey = ''
                    symm = False
                    continue
                if symm:
                    symmlist.append(line)
                if line[0] != "_":
                    loop_body = True
                # Loop header started, collecting keywords from head:
                if line[0] == "_" and atkey:
                    loophead_list.append(line)
                    continue
                # We are in a loop and the header ended, so we collect data:
                if loop_body and atkey:
                    loopitem = {}  # a line from the loop body, e.g. an atom
                    loop_data_line = delimit_line(line)
                    if cont:  # a continuation line
                        cont = False
                        continue
                    # unwrap loop data:
                    if len(loop_data_line) != len(loophead_list):
                        if textlen - 1 > num:
                            loop_data_line.extend(delimit_line(txt[num + 1].strip("\r\n ")))
                            cont = True
                        continue
                    for n, item in enumerate(loop_data_line):
                        loopitem[loophead_list[n]] = item
                    if cont:
                        continue
                    # TODO: make this general. Not only for atoms:
                    if atkey and loopitem[atkey] in atoms:
                        # atom is already there, upating values
                        atoms[loopitem[atkey]].update(loopitem)
                    elif atkey:
                        # atom is not there, creating key
                        atoms[loopitem[atkey]] = loopitem  # loopitem[atkey] is the atoms name
                continue
            # Leave out save_ frames:
            if save_frame:
                continue
            if line[:5] == "save_":
                save_frame = True
                continue
            elif line[:5] == "save_" and save_frame:
                save_frame = False
            # First find the start of the cif (the data tag)
            if line[:5] == 'data_':
                if line == "data_global":
                    continue
                if not data:
                    name = line.split('_')[1].strip('\n\r')
                    self.cif_data['data'] = name
                    data = True
                    continue
                else:  # break in occurence of a second data_
                    break
            # Find the loop positions:
            if line[:5] == "loop_":
                loop = True
                continue
            # Collect all data items outside loops:
            if line.startswith('_') and not loop:
                lsplit = line.split()
                # add regular cif items:
                if len(lsplit) > 1:
                    self.cif_data[lsplit[0]] = " ".join(delimit_line(" ".join(lsplit[1:])))
                    continue
                # add one-liners that are just in the next line:
                if len(lsplit) == 1 and txt[num + 1]:
                    if txt[num + 1][0] != ';' and txt[num + 1][0] != "_":
                        self.cif_data[lsplit[0]] = " ".join(delimit_line(txt[num + 1]))
                        continue
            if line.startswith("_shelx_hkl_file") or line.startswith("_refln_"):
                hkl = True
                continue
            # Leave out hkl frames:
            if hkl:
                break
                # continue  # use continue if data is behind hkl
            if line.lstrip()[0] == ";" and hkl:
                hkl = False
                continue
            if semi_colon_text_field:
                if not line.lstrip().startswith(";"):
                    semi_colon_text_list.append(line)
                if (textlen - 1 > num) and txt[num + 1][0] == ";":
                    self.cif_data[semi_colon_text_field] = "{}".format(os.linesep).join(semi_colon_text_list)
                    semi_colon_text_list.clear()
                    semi_colon_text_field = ''
                    continue
            if (textlen - 1 > num) and txt[num + 1][0] == ";":
                # if line.startswith("_shelx_res_file"):
                #    break
                # continue  # use continue if data is behind res file
                semi_colon_text_field = line
                continue
        self.cif_data['_atom'] = atoms
        self.cif_data['_space_group_symop_operation_xyz'] = '\n'.join(symmlist)
        self.cif_data['file_length_lines'] = num + 1
        self.symmcards = SymmCards()
        symmcards = [x.replace("'", "").replace(" ", "").split(',') for x in symmlist]
        for s in symmcards:
            self.symmcards.append(s)
        tmp = None
        for sym in self.symmcards:
            if sym == SymmetryElement(['X', 'Y', 'Z']):
                self.cif_data["_space_group_centring_type"] = 'P'
                continue
            if sym == SymmetryElement(['0.5', '0.5', '0.5']):
                self.cif_data["_space_group_centring_type"] = 'I'
                continue
            if sym == SymmetryElement(['1/3', '2/3', '2/3']):
                if not tmp:
                    tmp = 'R'
                continue
            if sym == SymmetryElement(['2/3', '1/3', '1/3']) and tmp == 'R':
                self.cif_data["_space_group_centring_type"] = 'R'
                continue
            if sym == SymmetryElement(['0.0', '0.5', '0.5']):  # F-centered
                if not tmp:
                    tmp = 'F'
                continue
            if sym == SymmetryElement(['0.5', '0.0', '0.5']):
                if not tmp:
                    tmp = 'F'
                continue
            if sym == SymmetryElement(['0.5', '0.5', '0.0']):
                if tmp == 'F':
                    self.cif_data["_space_group_centring_type"] = 'F'
                continue
            if sym == SymmetryElement(['0.0', '0.5', '0.5']):
                self.cif_data["_space_group_centring_type"] = 'A'
                continue
            if sym == SymmetryElement(['0.5', '0.0', '0.5']):
                self.cif_data["_space_group_centring_type"] = 'B'
                continue
            if sym == SymmetryElement(['0.5', '0.5', '0.0']):
                self.cif_data["_space_group_centring_type"] = 'C'
                continue
        if not data:
            return False
        # if not atoms:
        #    self.cif_data.clear()
        #    return False
        else:
            self.handle_deprecates()
            return True

    def handle_deprecates(self):
        """
        Makes the old and new cif values equal.
        """
        if "_symmetry_space_group_name_H-M" in self.cif_data:
            self.cif_data["_space_group_name_H-M_alt"] = self.cif_data["_symmetry_space_group_name_H-M"]
        if "_diffrn_measurement_device" in self.cif_data:
            self.cif_data["_diffrn_measurement_device_type"] = self.cif_data["_diffrn_measurement_device"]
        if "_refine_ls_shift/esd_max" in self.cif_data:
            self.cif_data["_refine_ls_shift/su_max"] = self.cif_data["_refine_ls_shift/esd_max"]
        if "_diffrn_measurement_device" in self.cif_data:
            self.cif_data["_diffrn_measurement_device_type"] = self.cif_data["_diffrn_measurement_device"]
        if '_symmetry_space_group_name_Hall' in self.cif_data:
            self.cif_data['_space_group_name_Hall'] = self.cif_data['_symmetry_space_group_name_Hall']
        if '_symmetry_Int_Tables_number' in self.cif_data:
            self.cif_data['_space_group_IT_number'] = self.cif_data['_symmetry_Int_Tables_number']
        if '_diffrn_reflns_av_sigmaI/netI' in self.cif_data:
            self.cif_data['_diffrn_reflns_av_unetI/netI'] = self.cif_data['_diffrn_reflns_av_sigmaI/netI']

    def __iter__(self) -> dict:
        """
        An iterable for the Cif object
        :return: cif entries
        """
        if self.ok:
            yield self.cif_data
        else:
            yield {}

    def __hash__(self):
        return hash(self.cif_data)

    def __getattr__(self, item, item_alt=''):
        """
        Returns an attribute of the cif data dictionary.
        """
        if item in self.cif_data or item_alt in self.cif_data:
            try:
                return self.cif_data[item]
            except KeyError:
                return self.cif_data[item_alt]
        else:
            return ''

    def __str__(self) -> str:
        """
        The string representation for print(self)
        """
        out = ''
        for item in self.cif_data:
            if item == '_atom':
                out += "Atoms:         \t\t\t"+str(len(self.cif_data['_atom']))+'\n'
                continue
            out += item+':  \t'+"'"+str(self.cif_data[item])+"'"+'\n'
        return out

    @property
    def cell(self):
        a = self.cif_data['_cell_length_a']
        b = self.cif_data['_cell_length_b']
        c = self.cif_data['_cell_length_c']
        alpha = self.cif_data['_cell_angle_alpha']
        beta = self.cif_data['_cell_angle_beta']
        gamma = self.cif_data['_cell_angle_gamma']
        if not all((a, b, c, alpha, beta, gamma)):
            return False
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
        return [a, b, c, alpha, beta, gamma]


def delimit_line(line: str) -> list:
    """
    Searches for delimiters in a cif line and returns a list of the respective values.
    >>> line = " 'C'  'C'   0.0033   0.0016   'some text inside' \\"more text\\""
    >>> delimit_line(line)
    ['C', 'C', '0.0033', '0.0016', 'some text inside', 'more text']
    >>> delimit_line("123  123 sdf")
    ['123', '123', 'sdf']
    >>> delimit_line("'-2  34' '234'")
    ['-2 34', '234']
    >>> delimit_line("'x, y, z'")
    ['x, y, z']
    >>> delimit_line("'a dog's life'")
    ["a dog's life"]
    """
    data = []
    line = line.split(' ')
    word = []
    cont = False
    for i in line:
        if i:
            if i[0] == "'" or i[0] == '"':
                if i[-1] == "'" or i[-1] == '"':
                    data.append(i.strip("'").strip('"'))
                    continue
            if i[0] == "'" or i[0] == '"':
                word.clear()
                word.append(i.strip("'").strip('"'))
                cont = True
                continue
            if i[-1] == "'" or i[-1] == '"':
                word.append(i.strip("'").strip('"'))
                data.append(' '.join(word))
                cont = False
                continue
            if cont:
                word.append(i)
            else:
                data.append(i)
    return data


if __name__ == '__main__':
    cif = Cif()
    with open('test-data/p-1_a.cif', 'r') as f:
        cifok = cif.parsefile(f.readlines())
    pprint(cif.cif_data)
    print(cifok)
    import doctest
    failed, attempted = doctest.testmod()  # verbose=True)
    if failed == 0:
        print('passed all {} tests!'.format(attempted))

