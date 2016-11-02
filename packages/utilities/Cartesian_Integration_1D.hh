#ifndef Cartesian_Integration_1D_hh
#define Cartesian_Integration_1D_hh

#include <vector>

#include "Integration_1D.hh"

class Cartesian_Integration_1D : public Integration_1D
{
public:
    
    Cartesian_Integration_1D(int order);
    
    double integrate(Integrand_1D &func,
                     double x1,
                     double x2) const;
    
private:
    
    int order_;
    std::vector<double> ordinates_;
    std::vector<double> weights_;
};

#endif
