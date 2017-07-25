#include "Reaction_Rates.hh"

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"

using namespace std;

Reaction_Rates::
Reaction_Rates(shared_ptr<Weak_Spatial_Discretization> spatial,
               shared_ptr<Angular_Discretization> angular,
               shared_ptr<Energy_Discretization> energy):
    Vector_Operator(),
    spatial_(spatial),
    angular_(angular),
    energy_(energy)
{
    int number_of_points = spatial->number_of_points();
    int number_of_moments = angular->number_of_moments();
    int number_of_groups = energy->number_of_groups();
    column_size_ = number_of_points * number_of_groups * number_of_moments;
    row_size_ = 3 * number_of_groups;
}

void Reaction_Rates::
apply(vector<double> &x) const
{
    vector<double> result(row_size_, 0);

    int number_of_points = spatial_->number_of_points();
    int number_of_moments = angular_->number_of_moments();
    int number_of_groups = energy_->number_of_groups();
    
    for (int i = 0; i < number_of_points; ++i)
    {
        // Get cross sections
        vector<double> sigma_t;
        vector<double> sigma_s;
        vector<double> sigma_f;

        get_coefficients(i,
                         sigma_t,
                         sigma_s,
                         sigma_f);

        for (int g = 0; g < number_of_groups; ++g)
        {
            // Add to reaction rate total
        }
    }
}

void Reaction_Rates::
get_coefficients(int i,
                 vector<double> &sigma_t,
                 vector<double> &sigma_s,
                 vector<double> &sigma_f) const
{
    int number_of_moments = angular_->number_of_moments();
    int number_of_scattering = angular_->number_of_scattering_moments();
    int number_of_groups = energy_->number_of_groups();
    bool identical_basis_functions
        = (spatial_->options()->identical_basis_functions
           == Weak_Spatial_Discretization_Options::Identical_Basis_Functions::TRUE);

    sigma_t.resize(number_of_groups);
    sigma_s.resize(number_of_groups);
    sigma_f.resize(number_of_groups);
    if (identical_basis_functions)
    {
        // Use weight function cross section integrals
        
    }
    else
    {
        // Use basis/weight function integrals with cross section integrals
        
    }
}
