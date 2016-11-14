#include "Cartesian_Integration_3D.hh"

#include "Quadrature_Rule.hh"

using std::vector;

Cartesian_Integration_3D::
Cartesian_Integration_3D(int order):
    order_(order)
{
    Quadrature_Rule::gauss_legendre(order,
                                    ordinates_,
                                    weights_);
}

double Cartesian_Integration_3D::
integrate(Integrand_3D &func,
          double x1,
          double x2,
          double y1,
          double y2,
          double z1,
          double z2) const
{
    double sum = 0;
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;
    double xt = x2 + x1;
    double yt = y2 + y1;
    double zt = z2 + z1;
    
    for (int i = 0; i < order_; ++i)
    {
        double point_x = 0.5 * (xt + dx * ordinates_[i]);

        for (int j = 0; j < order_; ++j)
        {
            double point_y = 0.5 * (yt + dy * ordinates_[j]);
            
            for (int k = 0; k < order_; ++k)
            {
                double point_z = 0.5 * (zt + dz * ordinates_[k]);
                
                sum += weights_[i] * weights_[j] * weights_[k] * func(point_x, point_y, point_z);
            }
        }
    }
    
    return 0.125 * dx * dy * dz * sum;
}
