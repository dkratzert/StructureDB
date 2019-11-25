/*cppimport
<%
setup_pybind11(cfg)
%>
*/
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
//#include <include/G6.h>
#include <include/NCDist.h>
//#include <include/LatticeConverter.h>
//#include <include/LRL_Cell.h>

namespace py = pybind11;


double NCDist_wrapper(std::vector<double> lst1, std::vector<double> lst2){
    // Wrapper method for python->C++ conversion
    /*
    5.2601, 9.1644, 10.609, 104.851, 104.324, 100.457
    5.2684, 9.208, 10.6641, 69.559, 76.132, 79.767
    raw: 49.6829
    0.7048613134033543
    */
    double mm1[6];
    double mm2[6];
    for (int i=0; i<6; i++) {
        mm1[i] = lst1[i];
        std::cout << mm1[i] << ", ";
    }
    std::cout << "\n";
    for (int j=0; j<6; j++) {
        mm2[j] = lst2[j];
        std::cout << mm2[j] << ", ";
    }
    std::cout << "\n";
    double raw_ncdist = NCDist(mm1, mm2);
    std::cout << "raw: " << raw_ncdist << "\n";
    return 0.1*std::sqrt(raw_ncdist);
}


/*G6 niggli_reduce_wrapper(const std::string &lattice, std::vector<double> &cell){
    LRL_Cell lcell = LRL_Cell(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]);
    G6 prim1 = LatticeConverter::NiggliReduceCell(lattice, lcell);
    return prim1;
}*/

PYBIND11_MODULE(ncdist, m) {
    m.def("ncdist", &NCDist_wrapper);
    //m.def("niggli_reduce", &niggli_reduce_wrapper);
}

