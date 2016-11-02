#ifndef Cylindrical_Integration_2D_hh
#define Cylindrical_Integration_2D_hh

#include <vector>

#include "Integration_2D.hh"

class Cylindrical_Integration_2D : public Integration_2D
{
public:

    Cylindrical_Integration_2D(int order);
    
    double integrate(Integrand_2D &func,
                     double x0,
                     double y0,
                     double rmax) const;
    
private:
    
    int order_;
    std::vector<double> ordinates_;
    std::vector<double> weights_;
};

#endif
