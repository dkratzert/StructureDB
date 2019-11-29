/*cppimport

*/
//#include <omp.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <NCDist.h>


namespace py = pybind11;


double NCDist_wrapper(std::vector<double> lst1, std::vector<double> lst2){
    // Wrapper method for python->C++ conversion
    double mm1[6];
    double mm2[6];
    for (int i=0; i<6; i++) {
        mm1[i] = lst1[i];
    }
    for (int j=0; j<6; j++) {
        mm2[j] = lst2[j];
    }
    double raw_ncdist = NCDist(mm1, mm2);
    return raw_ncdist;
}


PYBIND11_MODULE(ncdist, m) {
    m.def("ncdist", &NCDist_wrapper);
}

