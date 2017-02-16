#include "Meshless_Function.hh"
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
