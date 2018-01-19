#include "Wendland1_RBF.hh"

#include <cmath>
#include <string>

#include "Check.hh"

using std::pow;
using std::to_string;

Wendland1_RBF::
Wendland1_RBF(int order):
    order_(order)
{
}

double Wendland1_RBF::
value(double r) const
{
    r = std::abs(r);

    if (r <= 1)
    {
        switch(order_)
        {
        case 0:
            return 1 - r;
        case 1:
            return pow(1 - r, 3) * (1 + 3 * r);
        case 2:
            return pow(1 - r, 5) * (1 + 5 * r + 8 * r * r);
        case 3:
            return pow(1 - r, 7) * (1 + 7 * r + 19 * r * r + 21 * r * r * r);
        default:
            AssertMsg(false, "order \"" + to_string(order_) + "\" not implemented");
        }
    }
    return 0;
}

double Wendland1_RBF::
d_value(double r) const
{
    r = std::abs(r);

    if (r <= 1)
    {
        switch(order_)
        {
        case 0: // not smooth at r = 1
            return -1;
        case 1:
            return -12 * pow(-1 + r, 2) * r;
        case 2:
            return -14 * pow(-1 + r, 4) * r * (1 + 4 * r);
        case 3:
            return -6 * pow(-1 + r, 6) * r * (3 + 18 * r + 35 * r * r);
        default:
            AssertMsg(false, "order \"" + to_string(order_) + "\" not implemented");
        }
    }
    
    return 0;
}

double Wendland1_RBF::
dd_value(double r) const
{
    r = std::abs(r);
    
    if (r <= 1)
    {
        switch(order_)
        {
        case 0: 
            return 0;
        case 1:
            return -12 * (1 - 4 * r + 3 * r * r);
        case 2:
            return -14 * pow(-1 + r, 3) * (-1 - 3 * r + 24 * r * r);
        case 3:
            return -18 * pow(-1 + r, 5) * (-1 - 5 * r + 13 * r * r + 105 * r * r * r);
        default:
            AssertMsg(false, "order \"" + to_string(order_) + "\" not implemented");
        }
    }
    
    return 0;
}
