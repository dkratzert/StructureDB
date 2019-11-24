/*cppimport
<%
setup_pybind11(cfg)
%>
*/
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <include/NCDist.h>

namespace py = pybind11;


double NCDist_wrapper2(std::vector<double> lst1, std::vector<double> lst2){
    // Wrapper method for python->C++ conversion
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

PYBIND11_MODULE(ncdist, m) {
    m.def("ncdist", &NCDist_wrapper);
    m.def("niggli_reduce", &niggli_reduce);
}
