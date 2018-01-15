#include "Manufactured_Constant_Cross_Sections.hh"

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"

using namespace std;

Manufactured_Constant_Cross_Sections::
Manufactured_Constant_Cross_Sections(shared_ptr<Angular_Discretization> angular,
                                     shared_ptr<Energy_Discretization> energy,
                                     vector<double> const &sigma_t,
                                     vector<double> const &sigma_s):
    Manufactured_Cross_Sections(angular,
                                energy),
    sigma_t_(sigma_t),
    sigma_s_(sigma_s)
{
    int number_of_groups = energy_->number_of_groups();
    int number_of_scattering_moments = angular_->number_of_scattering_moments();
    Assert(sigma_t_.size() == number_of_groups);
    Assert(sigma_s_.size() == number_of_groups * number_of_groups * number_of_scattering_moments);
}

void Manufactured_Constant_Cross_Sections::
get_cross_sections(std::vector<double> const &position,
                   std::vector<double> &sigma_t,
                   std::vector<double> &sigma_s) const
{
    sigma_t = sigma_t_;
    sigma_s = sigma_s_;
}
