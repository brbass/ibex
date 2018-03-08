#include "Manufactured_Gaussian_Pincell.hh"

#include <cmath>

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Factory.hh"
#include "Cartesian_Distance.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"

using namespace std;

Manufactured_Gaussian_Pincell::
Manufactured_Gaussian_Pincell(shared_ptr<Angular_Discretization> angular,
                              shared_ptr<Energy_Discretization> energy,
                              vector<double> solution,
                              vector<double> aval,
                              vector<double> bval,
                              vector<double> cval,
                              vector<double> dval):
    Manufactured_Solution(angular,
                          energy),
    solution_(solution),
    aval_(aval),
    bval_(bval),
    cval_(cval),
    dval_(dval),
    distance_(make_shared<Cartesian_Distance>(angular->dimension()))
{
    int dimension = angular_->dimension();
    int number_of_moments = angular_->number_of_moments();
    int number_of_groups = energy_->number_of_groups();
    int solution_size = number_of_moments * number_of_groups;
    
    Assert(solution_.size() == solution_size);
    Assert(aval_.size() == solution_size);
    Assert(bval_.size() == solution_size);
    Assert(cval_.size() == solution_size);
    Assert(dval_.size() == solution_size);

    origin_.assign(dimension, 0.0);
}

vector<double> Manufactured_Gaussian_Pincell::
get_solution(vector<double> const &position) const
{
    int dimension = angular_->dimension();
    int number_of_moments = angular_->number_of_moments();
    int number_of_groups = energy_->number_of_groups();

    double dist = distance_->distance(position,
                                      origin_);
    
    vector<double> vals(number_of_moments * number_of_groups, 0);
    for (int m = 0; m < number_of_moments; ++m)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k = g + number_of_groups * m;

            double mult
                = (exp(aval_[k] * dist)
                   * (1 + bval_[k] * exp(cval_[k] * pow(dist + dval_[k], 2))));
            
            vals[k] = solution_[k] * mult;
        }
    }
    
    return vals;
}

vector<double> Manufactured_Gaussian_Pincell::
get_grad_solution(vector<double> const &position) const
{
    int dimension = angular_->dimension();
    int number_of_moments = angular_->number_of_moments();
    int number_of_groups = energy_->number_of_groups();

    double dist = distance_->distance(position,
                                      origin_);
    vector<double> grad_dist = distance_->gradient_distance(position,
                                                            origin_);
    vector<double> vals(dimension * number_of_moments * number_of_groups, 0);
    for (int d = 0; d < dimension; ++d)
    {
        for (int m = 0; m < number_of_moments; ++m)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k = g + number_of_groups * m;
                int kd = d + dimension * (g + number_of_groups * m);

                double adist = aval_[k] * dist;
                double ddist = dval_[k] + dist;
                double cdpow = cval_[k] * pow(ddist, 2);
                double mult
                    = (aval_[k] * exp(adist)
                       * (1 + bval_[k] * exp(cdpow) * grad_dist[d])
                       + 2 * bval_[k] * cval_[k] * exp(adist + cdpow)
                       * ddist * grad_dist[d]);
                
                vals[kd] = solution_[k] * mult;
            }
        }
    }
    
    return vals;
}
