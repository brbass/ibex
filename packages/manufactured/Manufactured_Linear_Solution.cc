#include "Manufactured_Linear_Solution.hh"

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Factory.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"

using namespace std;

Manufactured_Linear_Solution::
Manufactured_Linear_Solution(shared_ptr<Angular_Discretization> angular,
                             shared_ptr<Energy_Discretization> energy,
                             vector<double> origin,
                             vector<double> slope,
                             vector<double> solution):
    Manufactured_Solution(angular,
                          energy),
    origin_(origin),
    slope_(slope),
    solution_(solution)
{
    int dimension = angular_->dimension();
    int number_of_moments = angular_->number_of_moments();
    int number_of_groups = energy_->number_of_groups();
    Assert(solution_.size() == number_of_groups * number_of_moments);
    Assert(origin_.size() == dimension);
    Assert(slope_.size() == dimension);
}

vector<double> Manufactured_Linear_Solution::
get_solution(vector<double> const &position) const
{
    int dimension = angular_->dimension();
    int number_of_moments = angular_->number_of_moments();
    int number_of_groups = energy_->number_of_groups();

    double mult = 1;
    for (int d = 0; d < dimension; ++d)
    {
        mult += slope_[d] * (position[d] - origin_[d]);
    }
    
    vector<double> vals(number_of_moments * number_of_groups, 0);
    for (int m = 0; m < number_of_moments; ++m)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k = g + number_of_groups * m;
            vals[k] = solution_[k] * mult;
        }
    }
    
    return vals;
}

vector<double> Manufactured_Linear_Solution::
get_grad_solution(vector<double> const &position) const
{
    int dimension = angular_->dimension();
    int number_of_moments = angular_->number_of_moments();
    int number_of_groups = energy_->number_of_groups();

    vector<double> vals(dimension * number_of_moments * number_of_groups, 0);
    for (int d = 0; d < dimension; ++d)
    {
        double mult = slope_[d];
        for (int m = 0; m < number_of_moments; ++m)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k = g + number_of_groups * m;
                int kd = d + dimension * (g + number_of_groups * m);

                vals[kd] = solution_[k] * mult;
            }
        }
    }
    
    return vals;
}
