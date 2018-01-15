#include "Manufactured_Constant_Solution.hh"

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Factory.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"

using namespace std;

Manufactured_Constant_Solution::
Manufactured_Constant_Solution(shared_ptr<Angular_Discretization> angular,
                      shared_ptr<Energy_Discretization> energy,
                      vector<double> solution):
    Manufactured_Solution(angular,
                          energy),
    solution_(solution)
{
    int dimension = angular_->dimension();
    int number_of_moments = angular_->number_of_moments();
    int number_of_groups = energy_->number_of_groups();
    Assert(solution_.size() == number_of_groups * number_of_moments);
    grad_solution_.assign(number_of_groups * number_of_moments * dimension, 0);
}

vector<double> Manufactured_Constant_Solution::
get_solution(vector<double> const &position) const
{
    return solution_;
}

vector<double> Manufactured_Constant_Solution::
get_grad_solution(vector<double> const &position) const
{
    return grad_solution_;
}
