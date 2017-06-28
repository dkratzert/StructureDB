"""
MOl V3000 format
"""
import os

import elements
from searcher import misc
from searcher.database_handler import StructureTable


class MolFile():
    """
    This mol file writer is only to use the file with JSmol, not to implement the standard exactly!
    """
    def __init__(self, id: str, db: StructureTable, cell: tuple):
        self.db = db
        self.atoms = self.db.get_atoms_table(id, cell, cartesian=True)
        self.bonds = self.get_conntable_from_atoms()
        self.bondscount = len(self.bonds)
        self.atomscount = len(self.atoms)

    def header(self) -> str:
        """
        For JSmol, I don't need a facy header.
        """
        return "{}{}{}".format(os.linesep, os.linesep, os.linesep)

    def connection_table(self) -> str:
        """
          6  6  0  0  0  0  0  0  0  0  1 V3000
        """
        tab = "{0}  {1}  {2}  {3}  {4}  {5}  {6}  {7}  {8}  {9}  {10} {11}".format(self.atomscount,
                                            self.bondscount, 0, 0, 0, 0, 0, 0, 0, 0, "0999", "V2000")
        return tab

    def get_atoms_string(self) -> str:
        """
        Returns a string with an atom in each line.
        """
        atoms = []
        num = 0
        zeros = "0  0  0  0  0  0  0  0  0  0  0  0"
        for num, at in enumerate(self.atoms):
            x = at[2]
            y = at[3]
            z = at[4]
            element = at[1]
            atoms.append("{:>10.4f}{:>10.4f}{:>10.4f} {:<2s}".format(x, y, z, element, zeros))
        return '\n'.join(atoms)

    def get_bonds_string(self) -> str:
        """
        This is not like the file standard. The original standard wants to have fixed format 3 digitss for the bonds.
        """
        blist = []
        for bo in self.bonds:
            blist.append("{:>4d}{:>4d}  1  0  0  0  0".format(bo[0], bo[1]))
        return '\n'.join(blist)

    def get_conntable_from_atoms(self, extra_param = 0.16):
        """
        returns a connectivity table from the atomic coordinates and the covalence
        radii of the atoms.
        TODO:
        - read FREE command from db to contro binding here.
        :param cart_coords: cartesian coordinates of the atoms
        :type cart_coords: list
        :param atom_types: Atomic elements
        :type atom_types: list of strings
        :param atom_names: atom name in the file like C1
        :type atom_names: list of strings
        """
        conlist = []
        for num1, at1 in enumerate(self.atoms, 1):
            name1 = at1[0]
            typ1 = at1[1]
            x1 = at1[2]
            y1 = at1[3]
            z1 = at1[4]
            for num2, at2 in enumerate(self.atoms, 1):
                name2 = at2[0]
                typ2 = at2[1]
                x2 = at2[2]
                y2 = at2[3]
                z2 = at2[4]
                if name1 == name2:
                    continue
                d = misc.distance(x1, y1, z1, x2, y2, z2)
                # a bond is defined with less than the sum of the covalence
                # radii plus the extra_param:
                ele1 = elements.ELEMENTS[typ1.capitalize()]
                ele2 = elements.ELEMENTS[typ2.capitalize()]
                if d <= (ele1.covrad + ele2.covrad) + extra_param and d > (ele1.covrad or ele2.covrad):
                    conlist.append([num1, num2])
                    #print(num1, num2, d)
                    if len(conlist) == 99:
                        return conlist
                    if [num2, num1] or [num1, num2] in conlist:
                        continue
        return conlist

    def footer(self) -> str:
        """
        """
        return "M  END{}$$$$".format(os.linesep)

    def make_mol(self):
        """
        """
        header = '\n\n'
        connection_table = self.connection_table()
        atoms = self.get_atoms_string()
        bonds = self.get_bonds_string()
        footer = self.footer()
        mol = "{0}{5}{1}{5}{2}{5}{3}{5}{4}".format(header,connection_table,atoms,bonds,footer, '\n')
        return mol