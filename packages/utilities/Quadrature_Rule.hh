#ifndef Quadrature_Rule_hh
#define Quadrature_Rule_hh

#include <vector>

namespace Quadrature_Rule
{
    /*
      Get vectors of Gauss_Legendre ordinates and weights
    */
    void gauss_legendre(int n, std::vector<double> &ordinates, std::vector<double> &weights);
}
#endif


