#include "Compact_Gaussian_RBF.hh"

#include <cmath>

Compact_Gaussian_RBF::
Compact_Gaussian_RBF(double radius):
    radius_(radius)
{
    double k = exp(-radius_ * radius_);
    k1_ = 1 / (1 - k);
    k2_ = k / (1 - k);
}

double Compact_Gaussian_RBF::
value(double r) const
{
    if (r < radius_)
    {
        return k1_ * exp(-r * r) - k2_;
    }
    else
    {
        return 0.;
    }
}

double Compact_Gaussian_RBF::
d_value(double r) const
{
    if (r < radius_)
    {
        return -2 * k1_ * r * exp(-r * r);
    }
    else
    {
        return 0.;
    }
}

double Compact_Gaussian_RBF::
dd_value(double r) const
{
    if (r < radius_)
    {
        return 2 * (2 * r * r - 1) * k1_ * exp(-r * r);
    }
    else
    {
        return 0.;
    }
}
