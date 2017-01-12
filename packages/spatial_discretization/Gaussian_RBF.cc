#include "Gaussian_RBF.hh"

#include <cmath>

Gaussian_RBF::
Gaussian_RBF():
    RBF()
{
}

double Gaussian_RBF::
basis(double r) const
{
    return exp(-r * r);
}

double Gaussian_RBF::
d_basis(double r) const
{
    return -2 * r * exp(-r * r);
}

double Gaussian_RBF::
dd_basis(double r) const
{
    return (-2 + 4 * r * r) * exp(- r * r);
}
