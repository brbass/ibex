#include "Legendre.hh"

#include "Check.hh"
#include "Math_Functions.hh"

using namespace std;

Legendre::
Legendre(int order,
         vector<double> limits):
    order_(order),
    limits_(limits)
{
    Assert(limits_.size() == 2);
    length_ = limits_[1] - limits_[0];
    invlength_ = 1. / length_;
}

double Legendre::
get_local_position(double r) const
{
    return (2 * r - limits_[1] - limits_[0]) * invlength_;
}

double Legendre::
get_global_position(double xi) const
{
    return 0.5 * ((1 - xi) * limits_[0] + (1 + xi) * limits_[1]);
}

double Legendre::
value(double r) const
{
    double xi = get_local_position(r);
    return Math_Functions::legendre_polynomial(order_,
                                               xi);
}

double Legendre::
d_value(double r) const
{
    double xi = get_local_position(r);
    double val = Math_Functions::d_legendre_polynomial(order_,
                                                       xi);
    return 2 * invlength_ * val;
}

double Legendre::
dd_value(double r) const
{
    double xi = get_local_position(r);
    double val = Math_Functions::dd_legendre_polynomial(order_,
                                                        xi);
    return 4 * invlength_ * invlength_ * val;
}
