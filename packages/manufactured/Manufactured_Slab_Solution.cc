#include "Manufactured_Slab_Solution.hh"

#include <limits>

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Factory.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"

using namespace std;

Manufactured_Slab_Solution::
Manufactured_Slab_Solution(shared_ptr<Angular_Discretization> angular,
                           shared_ptr<Energy_Discretization> energy,
                           vector<double> interface_positions,
                           vector<vector<double> > solution):
    Manufactured_Solution(angular,
                          energy),
    number_of_regions_(solution.size()),
    interface_positions_(interface_positions),
    solution_(solution)
{
    int dimension = angular_->dimension();
    int number_of_moments = angular_->number_of_moments();
    int number_of_groups = energy_->number_of_groups();
    for (int i = 0; i < number_of_regions_; ++i)
    {
        Assert(solution_[i].size() == number_of_groups * number_of_moments);
    }
    grad_solution_.assign(number_of_groups * number_of_moments * dimension, 0);
    interface_positions_.push_back(numeric_limits<double>::max());
    Assert(interface_positions_.size() == number_of_regions_);
    for (int i = 0; i < number_of_regions_ - 1; ++i)
    {
        Assert(interface_positions_[i] < interface_positions_[i + 1]);
    }
}

vector<double> Manufactured_Slab_Solution::
get_solution(vector<double> const &position) const
{
    // Find appropriate cross sections
    for (int i = 0; i < number_of_regions_; ++i)
    {
        if (position[0] < interface_positions_[i])
        {
            return solution_[i];
        }
    }
    
    AssertMsg(false, "interface positions set incorrectly");
    return vector<double>();
}

vector<double> Manufactured_Slab_Solution::
get_grad_solution(vector<double> const &position) const
{
    return grad_solution_;
}
