#include "Heat_Transfer_Solution.hh"

#include "Check.hh"
#include "Weak_Spatial_Discretization.hh"

using namespace std;

Heat_Transfer_Solution::
Heat_Transfer_Solution(shared_ptr<Weak_Spatial_Discretization> spatial,
                       vector<double> const &coefficients):
    spatial_(spatial),
    coefficients_(coefficients)
{
    Assert(spatial);
}

double Heat_Transfer_Solution::
solution(std::vector<double> const &position) const
{
    Assert(position.size() == spatial_->dimension());
    
    return spatial_->expansion_value(position,
                                     coefficients_);
}

