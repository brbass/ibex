#include "Meshless_Function.hh"

#include "Check.hh"
#include "Vector_Functions.hh"

namespace vf = Vector_Functions;
using std::vector;

Meshless_Function::
Meshless_Function()
{
}

bool Meshless_Function::
inside_radius(vector<double> const &r) const
{
    double dist2 = vf::magnitude_squared(vf::subtract(r, position()));
    
    return dist2 < radius() * radius();
}

void Meshless_Function::
values(vector<double> const &r,
       vector<int> &indices,
       vector<double> &values) const
{
    AssertMsg(false, "not implemented for this meshless function");
}

void Meshless_Function::
gradient_values(vector<double> const &r,
                vector<int> &indices,
                vector<double> &vals,
                vector<vector<double> > &grad_vals) const
{
    AssertMsg(false, "not implemented for this meshless function");
}

