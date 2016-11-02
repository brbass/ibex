#ifndef Spherical_Integration_3D_hh
#define Spherical_Integration_3D_hh

#include <vector>

#include "Integration_3D.hh"

class Spherical_Integration_3D : public Integration_3D
{
public:
    
    Spherical_Integration_3D(int order);
    
    double integrate(Integrand_3D &func,
                     double x0,
                     double y0,
                     double z0,
                     double rmax) const;
    
    
private:
    
    int order_;
    std::vector<double> ordinates_;
    std::vector<double> weights_;
};

#endif
