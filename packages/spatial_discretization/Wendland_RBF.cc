#include "Wendland_RBF.hh"

#include <cmath>
#include <string>

#include "Check.hh"

using std::abs;
using std::pow;
using std::to_string;

Wendland_RBF::
Wendland_RBF(int order):
    RBF(),
    order_(order)
{
}

double Wendland_RBF::
basis(double r) const
{
    r = abs(r);

    if (r <= 1)
    {
        switch(order_)
        {
        case 0:
            return pow(1 - r, 2);
        case 1:
            return pow(1 - r, 4) * (1 + 4 * r);
        case 2:
            return pow(1 - r, 6) * (3 + 18 * r + 35 * r * r);
        case 3:
            return pow(1 - r, 8) * (1 + 8 * r + 25 * r * r + 32 * r * r * r);
        default:
            AssertMsg(false, "order \"" + to_string(order_) + "\" not implemented");
        }
    }
    return 0;
}

double Wendland_RBF::
d_basis(double r) const
{
    r = abs(r);

    if (r <= 1)
    {
        switch(order_)
        {
        case 0:
            return 2 * (-1 + r);
        case 1:
            return 20 * pow(-1 + r, 3) * r;
        case 2:
            return 56 * pow(-1 + r, 5) * r * (1 + 5 * r);
        case 3:
            return 22 * pow(-1 + r, 7) * r * (1 + 7 * r + 16 * r * r);
        default:
            AssertMsg(false, "order \"" + to_string(order_) + "\" not implemented");
        }
    }
    
    return 0;
}
    
double Wendland_RBF::
dd_basis(double r) const
{
    r = abs(r);
    
    if (r <= 1)
    {
        switch(order_)
        {
        case 0: // not smooth at r = 1
            return 2;
        case 1:
            return 20 * pow(-1 + r, 2) * (-1 + 4 * r);
        case 2:
            return 56 * pow(-1 + r, 4) * (-1 - 4 * r + 35 * r * r);
        case 3:
            return 22 * pow(-1 + r, 6) * (-1 - 6 * r + 15 * r * r + 160 * r * r * r);
        default:
            AssertMsg(false, "order \"" + to_string(order_) + "\" not implemented");
        }
    }
    
    return 0;
}
