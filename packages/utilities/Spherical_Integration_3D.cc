#include "Spherical_Integration_3D.hh"

#include <cmath>

#include "Quadrature_Rule.hh"

using std::vector;

Spherical_Integration_3D::
Spherical_Integration_3D(int order):
    order_(order)
{
    Quadrature_Rule::gauss_legendre(order,
                                    ordinates_,
                                    weights_);
}

double Spherical_Integration_3D::
integrate(Integrand_3D &func,
          double x0,
          double y0,
          double z0,
          double rmax) const
{
    double sum = 0;
    
    for (int i = 0; i < order_; ++i)
    {
        double point_r = 0.5 * (1 + ordinates_[i]) * rmax;

        for (int j = 0; j < order_; ++j)
        {
            double point_t = (1 + ordinates_[j]) * M_PI;
            
            for (int k = 0; k < order_; ++k)
            {
                double point_f = 0.5 * (1 + ordinates_[k]) * M_PI;
                double point_x = x0 + point_r * cos(point_t) * sin(point_f);
                double point_y = y0 + point_r * sin(point_t) * sin(point_f);
                double point_z = z0 + point_r * cos(point_f);
                
                sum += weights_[i] * weights_[j] * weights_[k] * func(point_x, point_y, point_z) * point_r * point_r * sin(point_f);
            }
        }
    }
    
    return 0.25 * M_PI * M_PI * rmax * sum;
}
