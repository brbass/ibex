#include "RBF.hh"

#include "Check.hh"

RBF::
RBF()
{
}

double RBF::
d_value(double s,
        double ds) const
{
    return d_value(s) * ds;
}

double RBF::
dd_value(double s,
         double ds2,
         double dds) const
{
    return dd_value(s) * ds2 + d_value(s) * dds;
}
