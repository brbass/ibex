#include "Inverse_Multiquadric_RBF.hh"

#include <cmath>

Inverse_Multiquadric_RBF::
Inverse_Multiquadric_RBF()
{
}

double Inverse_Multiquadric_RBF::
value(double r) const
{
    return 1 / sqrt(1 + r * r);
}

double Inverse_Multiquadric_RBF::
d_value(double r) const
{
    return -r * pow(1 + r * r, -1.5);
}

double Inverse_Multiquadric_RBF::
dd_value(double r) const
{
    return (-1 + 2 * r * r) * pow(1 + r * r, -2.5);
}
