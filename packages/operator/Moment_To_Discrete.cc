#include "Moment_To_Discrete.hh"

#include "Check.hh"
#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Spatial_Discretization.hh"

using namespace std;

Moment_To_Discrete::
Moment_To_Discrete(shared_ptr<Spatial_Discretization> spatial_discretization,
                   shared_ptr<Angular_Discretization> angular_discretization,
                   shared_ptr<Energy_Discretization> energy_discretization):
    Vector_Operator(),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization)
{
    int phi_size = (spatial_discretization->number_of_points()
                    * spatial_discretization->number_of_nodes()
                    * energy_discretization->number_of_groups()
                    * angular_discretization->number_of_moments());
    int psi_size = (spatial_discretization->number_of_points()
                    * spatial_discretization->number_of_nodes()
                    * energy_discretization->number_of_groups()
                    * angular_discretization->number_of_ordinates());
    row_size_ = psi_size;
    column_size_ = phi_size;
    
    check_class_invariants();
}

void Moment_To_Discrete::
apply(vector<double> &x) const
{
    vector<double> y(x);
    
    x.resize(row_size());
    
    int number_of_points = spatial_discretization_->number_of_points();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    double angular_normalization = angular_discretization_->angular_normalization();
    vector<int> const scattering_indices = angular_discretization_->scattering_indices();
    vector<double> const weights = angular_discretization_->weights();
    vector<double> const ordinates = angular_discretization_->ordinates();

    for (int i = 0; i < number_of_points; ++i)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            for (int n = 0; n < number_of_nodes; ++n)
            {
                for (int o = 0; o < number_of_ordinates; ++o)
                {
                    double sum = 0;
                    
                    for (int m = 0; m < number_of_moments; ++m)
                    {
                        int k = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                        int l = scattering_indices[m];

                        double p = angular_discretization_->moment(m, o);
                        
                        sum += (2 * static_cast<double>(l) + 1) / angular_normalization * p * y[k];
                    }
                    
                    int k = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                    
                    x[k] = sum;
                }
            }
        }
    }
}

void Moment_To_Discrete::
check_class_invariants() const
{
    Assert(spatial_discretization_);
    Assert(angular_discretization_);
    Assert(energy_discretization_);
}
