#include "Multiquadric_RBF.hh"

#include <cmath>

Multiquadric_RBF::
Multiquadric_RBF()
{
}

double Multiquadric_RBF::
value(double r) const
{
    return sqrt(1 + r * r);
}

double Multiquadric_RBF::
d_value(double r) const
{
    return r / sqrt(1 + r * r);
}

double Multiquadric_RBF::
dd_value(double r) const
{
    return pow(1 + r * r, -1.5);
}

