#ifndef Cartesian_Integration_3D_hh
#define Cartesian_Integration_3D_hh

#include <vector>

#include "Integration_3D.hh"

class Cartesian_Integration_3D : public Integration_3D
{
public:
    
    Cartesian_Integration_3D(int order);
    
    double integrate(Integrand_3D &func,
                     double x1,
                     double x2,
                     double y1,
                     double y2,
                     double z1,
                     double z2) const;
    
private:
    
    int order_;
    std::vector<double> ordinates_;
    std::vector<double> weights_;
};

#endif
