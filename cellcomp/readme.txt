                          ncdist -- Niggli Cone Distance utilities
                                 With S6 and D7 Support
                        Herbert J. Bernstein and Lawrence C. Andrews
                                Version 1.1.0 07 Jun 2019

                                    yayahjb@gmail.com


ncdist is a package of two utilities built on the Niggli-Cone cell distance discussed in
[Andrews, L. C.,  Bernstein, H. J. (2014). "The geometry of Niggli reduction: BGAOL --
embedding Niggli reduction and analysis of boundaries". J. Appl. Cryst.
47(1), 346 -- 359.] which provides a sensitive and reliable metric for the distance between
crystallographic unit cells.  Subsequent work has shown that similar distances can be
computed more efficiently based on Selling reduction which is the basis for Delaunay
reduction [Andrews, L. C., Bernstein, H. J. and Sauter, N. K. (2019). "Selling reduction
versus Niggli reduction for crystallographic lattices". Acta Cryst. A75(1), 115 -- 120.]
[Andrews, L. C., Bernstein, H. J. and Sauter, N. K. (2019). "A space for lattice
representation and clustering". Acta Cryst. A75(3), 593 -- 599.  Four shell command line
programs, ncdist, d7dist, s6dist_app, and cs6dist_app are provided as well as
a shared libraries, rcpp_ncdist (librcpp_ncdist.so), rcpp_d7dist (librcpp_d7dist.so),
rcpp_s6dist (librcpp_s6dist.so),  and rcpp_cs6dist (librcpp_cs6dist.so),  for use from R.

The 1.1.0 release of 7 Jun 2019 has been refactored to use code from L. C. Andrews
Lattice Representation Library (https://github.com/duck10/LatticeRepLib.git) with
changes to adapt to the conventions of the ncdist code.

Work on the upgrade to version 1.1.0 has been supported in part by funding from the US
Department of Energy Offices of Biological and Environmental Research and of Basic Energy
Sciences (grant No. DE-AC02-98CH10886 ; grant No. E-SC0012704); U.S. National Institutes
of Health (grant No. P41RR012408; grant No. P41GM103473; grant No. P41GM111244; and grant
No. R01GM117126 ); Dectris, Ltd.

LICENSE

#####################################################################
#                                                                    #
# YOU MAY REDISTRIBUTE THE ncdist PACKAGE UNDER THE TERMS OF THE GPL #
#                                                                    #
# ALTERNATIVELY YOU MAY REDISTRIBUTE THE ncdist API UNDER THE TERMS  #
# OF THE LGPL                                                        #
#                                                                    #
######################################################################

########################### GPL NOTICES ##############################
#                                                                    #
# This program is free software; you can redistribute it and/or      #
# modify it under the terms of the GNU General Public License as     #
# published by the Free Software Foundation; either version 2 of     #
# (the License, or (at your option) any later version.               #
#                                                                    #
# This program is distributed in the hope that it will be useful,    #
# but WITHOUT ANY WARRANTY; without even the implied warranty of     #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the      #
# GNU General Public License for more details.                       #
#                                                                    #
# You should have received a copy of the GNU General Public License  #
# along with this program; if not, write to the Free Software        #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA           #
# 02111-1307  USA                                                    #
#                                                                    #
######################################################################

######################### LGPL NOTICES ###############################
#                                                                    #
# This library is free software; you can redistribute it and/or      #
# modify it under the terms of the GNU Lesser General Public         #
# License as published by the Free Software Foundation; either       #
# version 2.1 of the License, or (at your option) any later version. #
#                                                                    #
# This library is distributed in the hope that it will be useful,    #
# but WITHOUT ANY WARRANTY; without even the implied warranty of     #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  #
# Lesser General Public License for more details.                    #
#                                                                    #
# You should have received a copy of the GNU Lesser General Public   #
# License along with this library; if not, write to the Free         #
# Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,    #
# MA  02110-1301  USA                                                #
#                                                                    #
######################################################################
