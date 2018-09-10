"""
MOl V3000 format
"""
import os

from sqlalchemy.orm import Session

from lattice import lattice
from searcher import misc
from searcher.atoms import get_radius_from_element
from searcher.database_handler import StructureTable, Atoms, Structure
from searcher.unitcell import Lattice
from shelxfile.dsrmath import Array


class MolFile(object):
    """
    This mol file writer is only to use the file with JSmol, not to implement the standard exactly!
    """
    def __init__(self, id: str, session: Session, cell: list, grow=False):
        self.session = session
        self.atoms = []
        if grow:
            pass
        else:
            atoms = [(at.Name, at.element, at.x, at.y, at.z) for at in
                     session.query(Atoms).filter(Atoms.StructureId == id).all() if at.Name]
            a = lattice.A(cell).orthogonal_matrix
            for at in atoms:
                self.atoms.append(at[:2] + tuple((Array([at[2], at[3], at[4]]) * a).values))
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
        tab = "{:>5d}{:>5d}".format(self.atomscount, self.bondscount)
        return tab

    def get_atoms_string(self) -> str:
        """
        Returns a string with an atom in each line.
        """
        atoms = []
        for num, at in enumerate(self.atoms):
            atoms.append("{:>10.4f}{:>10.4f}{:>10.4f} {:<2s}".format(at[2], at[3], at[4], at[1]))
        return '\n'.join(atoms)

    def get_bonds_string(self) -> str:
        """
        This is not accodingly to the file standard!
        The standard wants to have fixed format 3 digits for the bonds.
        """
        blist = []
        for bo in self.bonds:
            # This is deviating from the standard:
            blist.append("{:>4d}{:>4d}  1  0  0  0  0".format(bo[0], bo[1]))
        return '\n'.join(blist)

    def get_conntable_from_atoms(self, extra_param=0.35):
        """
        returns a connectivity table from the atomic coordinates and the covalence
        radii of the atoms.
        # a bond is defined with less than the sum of the covalence
        # radii plus the extra_param:
        TODO:
        - read FREE command from db to control binding here.
        :param extra_param: additional distance to the covalence radius
        :type extra_param: float
        """
        #t1 = time.clock()
        conlist = []
        for num1, at1 in enumerate(self.atoms, 1):
            rad1 = get_radius_from_element(at1[1])
            for num2, at2 in enumerate(self.atoms, 1):
                if at1[0] == at2[0]: # name1 = name2
                    continue
                rad2 = get_radius_from_element(at2[1])
                d = misc.distance(at1[2], at1[3], at1[4], at2[2], at2[3], at2[4])
                if (rad1 + rad2) + extra_param >= d > (rad1 or rad2):
                    conlist.append([num1, num2])
                    #print(num1, num2, d)
                    if [num2, num1] in conlist:
                        continue
        #t2 = time.clock()
        #print(round(t2-t1, 4), 's')
        return conlist

    def footer(self) -> str:
        """
        """
        return "M  END{}$$$$".format(os.linesep)

    def make_mol(self):
        """
        Combines all above to a mol file.
        """
        header = '\n\n'
        connection_table = self.connection_table()
        atoms = self.get_atoms_string()
        bonds = self.get_bonds_string()
        footer = self.footer()
        mol = "{0}{5}{1}{5}{2}{5}{3}{5}{4}".format(header,connection_table,atoms,bonds,footer, '\n')
        return mol