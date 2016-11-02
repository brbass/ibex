#ifndef Cartesian_Integration_2D_hh
#define Cartesian_Integration_2D_hh

#include <vector>

#include "Integration_2D.hh"

class Cartesian_Integration_2D : public Integration_2D
{
public:
    
    Cartesian_Integration_2D(int order);
    
    double integrate(Integrand_2D &func,
                     double x1,
                     double x2,
                     double y1,
                     double y2) const;
    
private:
    
    int order_;
    std::vector<double> ordinates_;
    std::vector<double> weights_;
};

#endif
