#include "Cartesian_Integration_1D.hh"

#include "Quadrature_Rule.hh"

using std::vector;

Cartesian_Integration_1D::
Cartesian_Integration_1D(int order):
    order_(order)
{
    Quadrature_Rule::gauss_legendre(order,
                                    ordinates_,
                                    weights_);
}

double Cartesian_Integration_1D::
integrate(Integrand_1D &func,
          double x1,
          double x2) const
{
    double sum = 0;
    double dx = x2 - x1;
    double xt = x2 + x1;
    
    for (int i = 0; i < order_; ++i)
    {
        double point = 0.5 * (xt + dx * ordinates_[i]);
        
        sum += weights_[i] * func(point);
    }

    return 0.5 * dx * sum;
}
