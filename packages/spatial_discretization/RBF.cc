#include "RBF.hh"

#include "Check.hh"

RBF::
RBF()
{
}

double RBF::
d_basis(double s,
        double ds) const
{
    return d_basis(s) * ds;
}

double RBF::
dd_basis(double s,
         double ds2,
         double dds) const
{
    return dd_basis(s) * ds2 + d_basis(s) * dds;
}
