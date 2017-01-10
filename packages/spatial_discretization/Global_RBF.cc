#include "Global_RBF.hh"

#include <limits>

Global_RBF::
Global_RBF()
{
}

double Global_RBF::
radius()
{
    return std::numeric_limits<double>::max();
}
