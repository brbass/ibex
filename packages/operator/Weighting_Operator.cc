#include "Weighting_Operator.hh"

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"

using std::shared_ptr;

Weighting_Operator::
Weighting_Operator(shared_ptr<Weak_Spatial_Discretization> spatial,
                   shared_ptr<Angular_Discretization> angular,
                   shared_ptr<Energy_Discretization> energy,
                   Options options):
    spatial_(spatial),
    angular_(angular),
    energy_(energy),
    options_(options)
{
}
