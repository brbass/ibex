#include "Manufactured_Radial_Cross_Sections.hh"

#include <cmath>
#include <limits>

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"

using namespace std;

Manufactured_Radial_Cross_Sections::
Manufactured_Radial_Cross_Sections(shared_ptr<Angular_Discretization> angular,
                                   shared_ptr<Energy_Discretization> energy,
                                   vector<double> const &interface_positions,
                                   vector<vector<double> > const &sigma_t,
                                   vector<vector<double> > const &sigma_s):
    Manufactured_Cross_Sections(angular,
                                energy),
    number_of_regions_(sigma_t.size()),
    interface_positions_(interface_positions),
    sigma_t_(sigma_t),
    sigma_s_(sigma_s)
{
    int number_of_groups = energy_->number_of_groups();
    int number_of_scattering_moments = angular_->number_of_scattering_moments();
    Assert(sigma_s_.size() == number_of_regions_);
    Assert(interface_positions_.size() == number_of_regions_ - 1);
    for (int i = 0; i < number_of_regions_; ++i)
    {
        Assert(sigma_t_[i].size() == number_of_groups);
        Assert(sigma_s_[i].size() == number_of_groups * number_of_groups * number_of_scattering_moments);
    }

    for (int i = 0; i < number_of_regions_ - 2; ++i)
    {
        Assert(interface_positions_[i] < interface_positions_[i + 1]);
        Assert(interface_positions_[i] > 0);
    }
}

void Manufactured_Radial_Cross_Sections::
get_cross_sections(std::vector<double> const &position,
                   std::vector<double> &sigma_t,
                   std::vector<double> &sigma_s) const
{
    // Get distance from origin
    int dimension = angular_->dimension();
    double dist = 0;
    for (int d = 0; d < dimension; ++d)
    {
        dist += position[d] * position[d];
    }
    dist = sqrt(dist);
    
    // Find appropriate cross sections
    for (int i = 0; i < number_of_regions_ - 1; ++i)
    {
        if (dist < interface_positions_[i])
        {
            sigma_t = sigma_t_[i];
            sigma_s = sigma_s_[i];
            return;
        }
    }

    sigma_t = sigma_t_[number_of_regions_ - 1];
    sigma_s = sigma_s_[number_of_regions_ - 1];
}
