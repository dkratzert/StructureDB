from math import pi, cos, sqrt, radians
from typing import List, Union

import numpy as np
from numpy.core.multiarray import ndarray
from spglib import spglib

# noinspection PyUnresolvedReferences
from cellcomp import ncdist
from shelxfile.misc import time_this_method


def abs_cap(val, max_abs_val=1):
    """
    Returns the value with its absolute value capped at max_abs_val.
    Particularly useful in passing values to trignometric functions where
    numerical errors may result in an argument > 1 being passed in.
    """
    return max(min(val, max_abs_val), -max_abs_val)


# Transformation matrix from the standardized unit cell to the primitive cell:
"""https://atztogo.github.io/spglib/definition.html#transformation-to-a-primitive-cell"""
Pa = np.array([[1., 0., 0.],
               [0., 0.5, -0.5],
               [0.0, 0.5, 0.5]], dtype=np.float64)
Pa.setflags(write=False)
Pc = np.array([[0.5, 0.5, 0],
               [-0.5, 0.5, 0],
               [0, 0, 1.0]], dtype=np.float64)
Pc.setflags(write=False)
Pr = np.array([[2.0 / 3.0, -1. / 3.0, -1.0 / 3.0],
               [1.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0],
               [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]], dtype=np.float64)
Pr.setflags(write=False)
Pi = np.array([[-0.5, 0.5, 0.5],
               [0.5, -0.5, 0.5],
               [0.5, 0.5, -0.5]], dtype=np.float64)
Pi.setflags(write=False)
Pf = np.array([[0.0, 0.5, 0.5],
               [0.5, 0.0, 0.5],
               [0.5, 0.5, 0.0]], dtype=np.float64)
Pf.setflags(write=False)

unitcell_type = List[float]
lattice_vect_type = List[Union[List[float], ndarray]]


class Lattice():
    def __init__(self, lattice_vect: lattice_vect_type, latt_type='P'):
        """
        A class to hold the lattice properties and calculate niggli cone distances to other lattices.
        :param lattice_vect: The lattice vectors [[a], [b], [c]]
        :param latt_type: Lattice type symbol of the current lattice
        """
        self.lattice_type = latt_type.upper()
        self.lattice_vect = np.array(lattice_vect, dtype=np.float64).reshape((3, 3))
        self.lattice_vect.setflags(write=False)

    def __repr__(self):
        return str(self.lattice_vect)

    def copy(self) -> 'Lattice':
        """Deep copy of self."""
        return Lattice(self.lattice_vect.copy(), self.lattice_type)

    @property
    def matrix(self) -> np.ndarray:
        """Copy of matrix representing the Lattice"""
        return self.lattice_vect

    @classmethod
    def from_parameters(cls, cell: List[Union[float, int]], latt_type: str = 'P') -> 'Lattice':
        """
        Create a Lattice using unit cell lengths and angles (in degrees).
        """
        a, b, c, alpha, beta, gamma = cell
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

    @property
    def primitive_lattice(self) -> 'Lattice':
        """
        Transforms the lattice to a primitive lattice.
        >>> from cellcomp.compare import Lattice
        >>> l = Lattice.from_parameters([1, 1, 1, 90,90,89], 'c')
        >>> l.primitive_lattice
        [[ 4.91273797e-01 -4.99923848e-01  0.00000000e+00]
         [ 5.08726203e-01  4.99923848e-01  6.12323400e-17]
         [ 0.00000000e+00  0.00000000e+00  1.00000000e+00]]
        """
        if self.lattice_type == 'A':
            return Lattice(np.dot(np.transpose(self.lattice_vect), Pa).T, 'P')
        elif self.lattice_type == 'C':
            return Lattice(np.dot(np.transpose(self.lattice_vect), Pc).T, 'P')
        elif self.lattice_type == 'R':
            return Lattice(np.dot(np.transpose(self.lattice_vect), Pr).T, 'P')
        elif self.lattice_type == 'H':
            return Lattice(np.dot(np.transpose(self.lattice_vect), Pr).T, 'P')
        elif self.lattice_type == 'I':
            return Lattice(np.dot(np.transpose(self.lattice_vect), Pi).T, 'P')
        elif self.lattice_type == 'F':
            return Lattice(np.dot(np.transpose(self.lattice_vect), Pf).T, 'P')
        else:
            return self.copy()

    def lengths(self) -> List[float]:
        """
        The a, b, c parameters of the unit cell.
        :return: A list with a, b, c
        """
        return np.sqrt(np.sum(self.lattice_vect ** 2, axis=1)).tolist()

    def angles(self) -> List[float]:
        """
        Returns the angles (alpha, beta, gamma) of the lattice.
        """
        angles = np.zeros(3)
        for i in range(3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            angles[i] = abs_cap(np.dot(self.lattice_vect[j], self.lattice_vect[k]) /
                                (self.lengths()[j] * self.lengths()[k]))
        angles = np.arccos(angles) * 180.0 / pi
        return angles.tolist()

    @property
    def unit_cell(self) -> unitcell_type:
        return self.lengths() + self.angles()

    @property
    def G6(self) -> List[float]:
        """ Take a reduced Niggli Cell, and turn it into the G6 representation """
        uc = self.unit_cell
        a = uc[0] ** 2
        b = uc[1] ** 2
        c = uc[2] ** 2
        d = 2 * uc[1] * uc[2] * cos(radians(uc[3]))
        e = 2 * uc[0] * uc[2] * cos(radians(uc[4]))
        f = 2 * uc[0] * uc[1] * cos(radians(uc[5]))
        return [a, b, c, d, e, f]

    # @time_this_method
    def ncdist_fromcell(self, cell2: List[Union[float, int]], unitcell_type: str) -> float:
        """
        Does a niggli reduction, G6 vector and distance calculation from two given unit cells.
        :param cell2: Second unit cell to compare self to
        :param unitcell_type: Lattice centering type of the unit cell, e.g. P, R, or I

        The square root of the BGAOL Niggli cone embedding distance NCDist based on
        [a2, b2, c2, 2bccos(α), 2accos(β), 2abcos(γ)] with the distances scaled by 1/√6 and divided by the 
        reciprocal of the average length of cell edges f. The square root linearizes the metric to Angstrom units.
        """
        # cell->lattice_vec->make_primitive->niggli_reduce->to_g6->ncdist_from_g6
        primlatt = Lattice.from_parameters(cell2, latt_type=unitcell_type).primitive_lattice
        dist = ncdist.ncdist(self.primitive_lattice.niggli_reduced.G6, primlatt.niggli_reduced.G6)
        # Is this really correct?
        dist = dist * (1 / sqrt(6))
        return 0.1 * sqrt(dist)

    def match_cell(self, cell2: List[Union[float, int]], unitcell_type: str) -> Union[int, float]:
        dist = self.ncdist_fromcell(cell2, unitcell_type)
        if dist < 0.2:
            return dist
        else:
            return 0


if __name__ == '__main__':
    c1 = [5.2601, 9.1644, 10.6090, 104.851, 104.324, 100.457]
    c2 = [5.2684, 9.2080, 10.6641, 69.559, 76.132, 79.767]
    c3 = [float(x) for x in "100 100 100 90 90 90".split()]
    c4 = [float(x) for x in "99 99 99 89 89 89".split()]
    c3a = [float(x) for x in "7.5675 13.1966 11.3486   90.000  103.608   90.000 ".split()]
    c4a = [float(x) for x in "7.6870 13.2020 11.5790   90.000  105.840   90.000 ".split()]

    c5 = [57.98, 57.98, 57.98, 92.02, 92.02, 92.02]  # R32  1FE5 3.3
    c6 = [80.36, 80.36, 99.44, 90.0, 90.0, 120.0]  # R3   1U4J, primitive= 57.02, 57.02, 57.02, 89.605, 89.605, 89.605
    c7 = [80.949, 80.572, 57.098, 90.0, 90.35, 90.0]  # C2  1G2X 0.9
    c6x = [77.516, 77.516, 99.076, 90.00, 90.00, 120.00]  # H3  2WCE 3
    c6y = [80.360, 80.360, 99.440, 90.00, 90.00, 120.00]  # H 1G0Z 0
    #
    c8 = [78.961, 82.328, 57.031, 90.00, 93.44, 90.00]  # C2  1GUT
    c9 = [80.36, 80.36, 99.44, 90, 90, 120]
    #
    c10 = [3.1457, 3.1457, 3.1451, 60.089, 60.0887, 60.104]
    c11 = [3.1456, 3.1458, 3.1451, 90.089, 119.907, 119.898]
    #
    # primitice cell:
    c12 = [10., 10., 10., 112., 112.5, 112.9]  # G6: (100, 100, 100, -74.9, -76.5, -77.8)
    # reduced cell:
    c13 = [8.41, 10.0, 10.0, 112, 106.3, 106.8]  # reduced G6: (70.7, 100, 100,-74.9,-47-3,-48.5)

    lat = Lattice.from_parameters(c6, 'R')
    dist = lat.ncdist_fromcell(c6y, 'h')
    print(dist, '1G0Z')
    dist = lat.ncdist_fromcell(c7, 'C')
    print(dist, '1G2X')
    dist = lat.ncdist_fromcell(c6x, 'h')
    print(dist, '2WCE')
    dist = lat.ncdist_fromcell(c8, 'c')
    print(dist, '1GUT')

    lat = Lattice.from_parameters(c1, 'p')
    dist = lat.ncdist_fromcell(c2, 'p')
    print(dist, 'wp')

    lat = Lattice.from_parameters(c3a, 'p')
    dist = lat.ncdist_fromcell(c4a, 'p')
    print(dist, 'wp2')

    lat = Lattice.from_parameters([18.7334, 25.9327, 7.3674, 90.000, 90.000, 90.000], 'p')
    dist = lat.ncdist_fromcell([7.3518 ,  18.6898,   25.8702 ,  90.0000  , 90.0000 ,  90.0000], 'p')
    print(dist, 'test')

    lat = Lattice.from_parameters([7.3518 ,  18.6898,   25.8702 ,  90.0000  , 90.0000 ,  90.0000], 'p')
    dist = lat.ncdist_fromcell([18.7334, 25.9327, 7.3674, 90.000, 90.000, 90.000], 'p')
    print(dist, 'test')

    @time_this_method
    def speed():
        for _ in range(1900):
            lat = Lattice.from_parameters(c3a, 'c')
            dist = lat.ncdist_fromcell(c4a, 'c')


    #speed()
