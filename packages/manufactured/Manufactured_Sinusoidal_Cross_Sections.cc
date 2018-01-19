#include "Manufactured_Sinusoidal_Cross_Sections.hh"

#include <cmath>
#include <iostream>

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"

using namespace std;

Manufactured_Sinusoidal_Cross_Sections::
Manufactured_Sinusoidal_Cross_Sections(shared_ptr<Angular_Discretization> angular,
                                       shared_ptr<Energy_Discretization> energy,
                                       double const relative_amplitude,
                                       vector<double> const &frequency,
                                       vector<double> const &sigma_t,
                                       vector<double> const &sigma_s):
    Manufactured_Cross_Sections(angular,
                                energy),
    relative_amplitude_(relative_amplitude),
    frequency_(frequency),
    sigma_t_(sigma_t),
    sigma_s_(sigma_s)
{
    int dimension = angular_->dimension();
    int number_of_groups = energy_->number_of_groups();
    int number_of_scattering_moments = angular_->number_of_scattering_moments();
    Assert(sigma_t_.size() == number_of_groups);
    Assert(sigma_s_.size() == number_of_groups * number_of_groups * number_of_scattering_moments);
    Assert(frequency_.size() == dimension);
    if (abs(relative_amplitude) > 1)
    {
        cout << "Manufactured_Sinusoidal_Cross_Sections: given relative amplitude may produce negative cross sections" << endl;
    }
}

void Manufactured_Sinusoidal_Cross_Sections::
get_cross_sections(std::vector<double> const &position,
                   std::vector<double> &sigma_t,
                   std::vector<double> &sigma_s) const
{
    int number_of_groups = energy_->number_of_groups();
    int number_of_scattering_moments = angular_->number_of_scattering_moments();
    int dimension = angular_->dimension();

    // Get multiplication factor
    double mult = relative_amplitude_;
    for (int d = 0; d < dimension; ++d)
    {
        mult *= cos(2 * M_PI * frequency_[d] * position[d]);
    }
    mult += 1;
    
    // Get new cross sections
    sigma_t = sigma_t_;
    sigma_s = sigma_s_;
    for (double &sigma : sigma_t)
    {
        sigma *= mult;
    }
    for (double &sigma : sigma_s)
    {
        sigma *= mult;
    }
}
