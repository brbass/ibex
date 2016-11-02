#include "Cylindrical_Integration_2D.hh"

#include <cmath>

#include "Quadrature_Rule.hh"

using std::vector;

Cylindrical_Integration_2D::
Cylindrical_Integration_2D(int order):
    order_(order)
{
    Quadrature_Rule::gauss_legendre(order,
                                    ordinates_,
                                    weights_);
}

double Cylindrical_Integration_2D::
integrate(Integrand_2D &func,
          double x0,
          double y0,
          double rmax) const
{
    double sum = 0;
    
    for (int i = 0; i < order_; ++i)
    {
        double point_r = 0.5 * (1 + ordinates_[i]) * rmax;
        
        for (int j = 0; j < order_; ++j)
        {
            double point_t = M_PI * (1 + ordinates_[j]);
            double point_x = x0 + point_r * cos(point_t);
            double point_y = y0 + point_r * sin(point_t);
            
            sum += weights_[i] * weights_[j] * func(point_x, point_y) * point_r;
        }
    }
    
    return 0.5 * M_PI * rmax * sum;
}
