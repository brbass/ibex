#include "Cartesian_Integration_2D.hh"

#include "Quadrature_Rule.hh"

using std::vector;

Cartesian_Integration_2D::
Cartesian_Integration_2D(int order):
    order_(order)
{
    Quadrature_Rule::gauss_legendre(order,
                                    ordinates_,
                                    weights_);
}

double Cartesian_Integration_2D::
integrate(Integrand_2D &func,
          double x1,
          double x2,
          double y1,
          double y2) const
{
    double sum = 0;
    double dx = x2 - x1;
    double dy = y2 - y1;
    double xt = x2 + x1;
    double yt = y2 + y1;
    
    for (int i = 0; i < order_; ++i)
    {
        double point_x = 0.5 * (xt + dx * ordinates_[i]);

        for (int j = 0; j < order_; ++j)
        {
            double point_y = 0.5 * (yt + dy * ordinates_[j]);
            
            sum += weights_[i] * weights_[j] * func(point_x, point_y);
        }
    }
    
    return 0.25 * dx * dy * sum;
}

