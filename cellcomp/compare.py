from math import pi, cos, sqrt, radians
from typing import List, Union, Iterable

import numpy as np
from spglib import spglib

# noinspection PyUnresolvedReferences
from cellcomp import ncdist


def abs_cap(val, max_abs_val=1):
    """
    Returns the value with its absolute value capped at max_abs_val.
    Particularly useful in passing values to trignometric functions where
    numerical errors may result in an argument > 1 being passed in.
    """
    return max(min(val, max_abs_val), -max_abs_val)


# Transformation matrix from the standardized unit cell to the primitive cell:
"""https://atztogo.github.io/spglib/definition.html#transformation-to-a-primitive-cell"""
Pa = np.array([[1, 0, 0],
               [0, 0.5, -0.5],
               [0, 0.5, 0.5]])
Pc = np.array([[0.5, 0.5, 0],
               [-0.5, 0.5, 0],
               [0, 0, 1.]])
Pr = np.array([[2. / 3., -1. / 3., -1 / 3.],
               [1 / 3., 1 / 3., -2 / 3.],
               [1 / 3., 1 / 3., 1 / 3.]])
Pi = np.array([[-0.5, 0.5, 0.5],
               [0.5, -0.5, 0.5],
               [0.5, 0.5, -0.5]])
Pf = np.array([[0, 0.5, 0.5],
               [0.5, 0, 0.5],
               [0.5, 0.5, 0]])

unitcell_type = Iterable[float]
lattice_vect_type = List[float]


class Lattice():
    def __init__(self, lattice_vect: lattice_vect_type, latt_type='P'):
        """
        A class to hold the lattice properties and calculate niggli cone distances to other lattices.
        :param lattice_vect: The lattice vectors [[a], [b], [c]]
        :param latt_type: Lattice type symbol of the current lattice
        """
        self.lattice_type = latt_type
        self.lattice_vect = np.array(lattice_vect, dtype=np.float64).setflags(write=False)

    @classmethod
    def from_parameters(cls, a: float, b: float, c: float, alpha: float, beta: float, gamma: float,
                        latt_type: str = 'P') -> 'Lattice':
        """
        Create a Lattice using unit cell lengths and angles (in degrees).
        """
        angles_r = np.radians([alpha, beta, gamma])
        cos_alpha, cos_beta, cos_gamma = np.cos(angles_r)
        sin_alpha, sin_beta, sin_gamma = np.sin(angles_r)
        val = abs_cap((cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta))
        gamma_star = np.arccos(val)
        vector_a = [a * sin_beta, 0.0, a * cos_beta]
        vector_b = [-b * sin_alpha * np.cos(gamma_star),
                    b * sin_alpha * np.sin(gamma_star),
                    b * cos_alpha, ]
        vector_c = [0.0, 0.0, float(c)]
        return Lattice([vector_a, vector_b, vector_c], latt_type)

    @property
    def niggli_reduced(self) -> 'Lattice':
        """
        Returns the niggli reduced basis vectors.
        """
        return Lattice(spglib.niggli_reduce(self.lattice_vect), self.lattice_type)

    def to_primitive_lattice(self) -> lattice_vect_type:
        """
        Transforms the lattice to a primitive lattice.
        """
        if self.lattice_type == 'A':
            return np.dot(np.transpose(self.lattice_vect), Pa).T
        elif self.lattice_type == 'C':
            return np.dot(np.transpose(self.lattice_vect), Pc).T
        elif self.lattice_type == 'R':
            return np.dot(np.transpose(self.lattice_vect), Pc).T
        elif self.lattice_type == 'I':
            return np.dot(np.transpose(self.lattice_vect), Pc).T
        elif self.lattice_type == 'F':
            return np.dot(np.transpose(self.lattice_vect), Pc).T
        else:
            return self.lattice_vect

    @staticmethod
    def lengths(matrix: lattice_vect_type) -> List[float]:
        """
        The a, b, c parameters of the unit cell.
        :param matrix: The lattice vector matrix.
        :return: A list with a, b, c
        """
        if not isinstance(matrix, np.ndarray):
            matrix = np.ndarray(matrix, dtype=float)
        return np.sqrt(np.sum(matrix ** 2, axis=1)).tolist()

    def angles(self, matrix: List[List[float]]) -> List[float]:
        """
        Returns the angles (alpha, beta, gamma) of the lattice.
        """
        lengt = self.lengths(matrix)
        angles = np.zeros(3)
        for i in range(3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            angles[i] = abs_cap(np.dot(matrix[j], matrix[k]) / (lengt[j] * lengt[k]))
        angles = np.arccos(angles) * 180.0 / pi
        return angles.tolist()

    def lattice_to_cell(self, lattice: lattice_vect_type) -> unitcell_type:
        return self.lengths(lattice) + self.angles(lattice)

    @staticmethod
    def cell_to_g6(uc: unitcell_type) -> unitcell_type:
        """ Take a reduced Niggli Cell, and turn it into the G6 representation """
        a = uc[0] ** 2
        b = uc[1] ** 2
        c = uc[2] ** 2
        d = 2 * uc[1] * uc[2] * cos(radians(uc[3]))
        e = 2 * uc[0] * uc[2] * cos(radians(uc[4]))
        f = 2 * uc[0] * uc[1] * cos(radians(uc[5]))
        return [a, b, c, d, e, f]

    def ncdist_fromcell(self, cell2: unitcell_type) -> float:
        """
        Does a niggli reduction, G6 vector and distance calculation from two given unit cells.
        TODO: transform to primitive lattices
        """
        reduced1 = self.niggli_reduced
        G6a = self.cell_to_g6(reduced1.to_primitive_lattice())
        reduced2 = Lattice.from_parameters(*cell2)
        G6b = self.cell_to_g6(reduced2.to_primitive_lattice())
        return 0.01 * sqrt(ncdist.ncdist(G6a, G6b))


if __name__ == '__main__':
    c1 = [5.2601, 9.1644, 10.6090, 104.851, 104.324, 100.457]
    c2 = [5.2684, 9.2080, 10.6641, 69.559, 76.132, 79.767]
    c3 = [float(x) for x in "100 100 100 90 90 90".split()]
    c4 = [float(x) for x in "99 99 99 89 89 89".split()]
    c3a = [float(x) for x in "7.5675 13.1966 11.3486   90.000  103.608   90.000 ".split()]
    c4a = [float(x) for x in "7.6870 13.2020 11.5790   90.000  105.840   90.000 ".split()]

    c5 = [57.98, 57.98, 57.98, 92.02, 92.02, 92.02]  # R32  1FE5
    c6 = [80.36, 80.36, 99.44, 90, 90, 120]  # R3   1U4J, primitive= 57.02, 57.02, 57.02, 89.605, 89.605, 89.605
    c7 = [80.949, 80.572, 57.098, 90.0, 90.35, 90.0]  # C2  1G2X
    #
    c8 = [78.961, 82.328, 57.031, 90.00, 93.44, 90.00]
    c9 = [80.36, 80.36, 99.44, 90, 90, 120]
    #
    c10 = [3.1457, 3.1457, 3.1451, 60.089, 60.0887, 60.104]
    c11 = [3.1456, 3.1458, 3.1451, 90.089, 119.907, 119.898]
    #
    # primitice cell:
    c12 = [10., 10., 10., 112., 112.5, 112.9]  # G6: (100, 100, 100, -74.9, -76.5, -77.8)
    # reduced cell:
    c13 = [8.41, 10.0, 10.0, 112, 106.3, 106.8]  # reduced G6: (70.7, 100, 100,-74.9,-47-3,-48.5)

    lat = Lattice(c5)
    dist = lat.ncdist_fromcell(c6)
    print(dist)
    """
    lat1 = Lattice(*c12)
    reduced = lattice_to_cell(spglib.niggli_reduce(lattice_from_parameters(*c12)))
    print('g6 of 12:', cell_to_g6(c12))
    print('reduced:', reduced)
    G6a = cell_to_g6(reduced)
    print('reduced g6 of 12:', G6a)

    print(ncdist_fromcell(c1, c2))
    print(ncdist_fromcell(c3a, c4a))
    print(ncdist_fromcell(c3, c4))
    print('c5:')
    print(ncdist_fromcell(c6, c5))
    print(ncdist_fromcell(c5, c7))
    # print(np.dot(np.transpose(lattice), Pc).T)
    reduced = spglib.niggli_reduce(lattice_from_parameters(*c6))
    print(lattice_to_cell(np.dot(np.transpose(reduced), Pr).T)) 
    """
