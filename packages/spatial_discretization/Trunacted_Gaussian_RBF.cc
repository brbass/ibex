#include "Truncated_Gaussian_RBF.hh"

#include <cmath>

Truncated_Gaussian_RBF::
Truncated_Gaussian_RBF(double radius):
    radius_(radius)
{
    double k = exp(-radius_ * radius_);
    k1_ = 1 / (1 - k);
    k2_ = k / (1 - k);
}

double Truncated_Gaussian_RBF::
basis(double r) const
{
    if (r < radius)
    {
        return exp(-r * r);
    }
    else
    {
        return 0.
    }
}

double Truncated_Gaussian_RBF::
d_basis(double r) const
{
    if (r < radius)
    {
        return -2 * r * exp(-r * r);
    }
    else
    {
        return 0.
    }
}

double Truncated_Gaussian_RBF::
dd_basis(double r) const
{
    if (r < radius)
    {
        return (-2 + 4 * r * r) * exp(- r * r);
    }
    else
    {
        return 0.
    }
}
