#include "Global_RBF.hh"

#include <limits>

Global_RBF::
Global_RBF()
{
}

double Global_RBF::
radius() const
{
    return 0.5 * std::numeric_limits<double>::max();
}
