/*cppimport

*/
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <NCDist.h>

namespace py = pybind11;

int square(int x) {
    return x * x;
}



double NCDist_wrapper(double a, double b, double c, double al, double be, double ga,
                        double a1, double b1, double c1, double al1, double be1, double ga1){
    double mm1[6];
    double mm2[6];
    mm1[0] = a;
    mm1[1] = b;
    mm1[2] = c;
    mm1[3] = al;
    mm1[4] = be;
    mm1[5] = ga;
    mm2[0] = a1;
    mm2[1] = b1;
    mm2[2] = c1;
    mm2[3] = al1;
    mm2[4] = be1;
    mm2[5] = ga1;
    return 0.1*std::sqrt(NCDist(mm1, mm2));
}

double NCDist_wrapper2(std::vector<double> lst1, std::vector<double> lst2){
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
    return 0.1*std::sqrt(NCDist(mm1, mm2));
}

PYBIND11_MODULE(ncdist, m) {
    m.def("ncdist", &NCDist_wrapper);
    m.def("ncdist2", &NCDist_wrapper2);
}
