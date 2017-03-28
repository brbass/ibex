#include "Linf_Convergence.hh"

#include <cmath>

#include "Check.hh"

using std::abs;
using std::vector;

Linf_Convergence::
Linf_Convergence(double tolerance):
    tolerance_(tolerance)
{
}

double Linf_Convergence::
error(double val_new,
      double val_old) const
{
    return abs(val_new - val_old);
}

double Linf_Convergence::
error(vector<double> const &val_new,
      vector<double> const &val_old) const
{
    int size = val_new.size();
    Assert(val_old.size() == size);
    
    double max = 0;
    for (int i = 0; i < size; ++i)
    {
        double err = abs(val_new[i] - val_old[i]);
        
        if (err > max)
        {
            max = err;
        }
    }

    return max;
}

bool Linf_Convergence::
check(double error_new,
      double error_old) const
{
    double rho = error_new / error_old;
    double ep = tolerance_ * (1 - rho);
    
    return error_new < ep;
}
