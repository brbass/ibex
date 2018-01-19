#include "Manufactured_Sinusoidal_Solution.hh"

#include <cmath>
#include <iomanip>
#include <iostream>

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"

using namespace std;

Manufactured_Sinusoidal_Solution::
Manufactured_Sinusoidal_Solution(shared_ptr<Angular_Discretization> angular,
                                 shared_ptr<Energy_Discretization> energy,
                                 double const relative_amplitude,
                                 vector<double> const &frequency,
                                 vector<double> const &solution):
    Manufactured_Solution(angular,
                          energy),
    relative_amplitude_(relative_amplitude),
    frequency_(frequency),
    solution_(solution)
{
    int dimension = angular_->dimension();
    int number_of_moments = angular_->number_of_moments();
    int number_of_groups = energy_->number_of_groups();
    Assert(solution_.size() == number_of_groups * number_of_moments);
    Assert(frequency_.size() == dimension);
    if (abs(relative_amplitude) > 1)
    {
        cout << "Manufactured_Sinusoidal_Solution: given relative amplitude may produce negative solution" << endl;
    }
}

vector<double> Manufactured_Sinusoidal_Solution::
get_solution(vector<double> const &position) const
{
    int number_of_groups = energy_->number_of_groups();
    int dimension = angular_->dimension();

    // Get multiplication factor
    double mult = relative_amplitude_;
    for (int d = 0; d < dimension; ++d)
    {
        mult *= cos(2 * M_PI * frequency_[d] * position[d]);
    }
    mult += 1;
    
    vector<double> vals = solution_;
    for (double &val : vals)
    {
        val *= mult;
    }
    
    return vals;
}

vector<double> Manufactured_Sinusoidal_Solution::
get_grad_solution(vector<double> const &position) const
{
    int number_of_groups = energy_->number_of_groups();
    int number_of_moments = angular_->number_of_moments();
    int dimension = angular_->dimension();

    // Get multiplication factor
    vector<double> mult(dimension, relative_amplitude_);
    for (int d1 = 0; d1 < dimension; ++d1)
    {
        for (int d2 = 0; d2 < dimension; ++d2)
        {
            if (d1 == d2)
            {
                mult[d1] *= -2 * M_PI * frequency_[d1] * sin(2 * M_PI * frequency_[d1] * position[d1]);
            }
            else
            {
                mult[d1] *= cos(2 * M_PI * frequency_[d2] * position[d2]);
            }
        }
    }
    
    vector<double> vals(number_of_groups * number_of_moments * dimension);
    for (int m = 0; m < number_of_moments; ++m)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            for (int d = 0; d < dimension; ++d)
            {
                int k = g + number_of_groups * m;
                int kd = d + dimension * (g + number_of_groups * m);
                
                vals[kd] = solution_[k] * mult[d];
            }
        }
    }
    
    return vals;
}
