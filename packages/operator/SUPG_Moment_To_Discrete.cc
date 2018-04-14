#include "SUPG_Moment_To_Discrete.hh"

#if defined(ENABLE_OPENMP)
    #include <omp.h>
#endif

#include "Check.hh"
#include "Angular_Discretization.hh"
#include "Dimensional_Moments.hh"
#include "Energy_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"

using namespace std;

SUPG_Moment_To_Discrete::
SUPG_Moment_To_Discrete(shared_ptr<Weak_Spatial_Discretization> spatial_discretization,
                        shared_ptr<Angular_Discretization> angular_discretization,
                        shared_ptr<Energy_Discretization> energy_discretization,
                        bool include_double_dimensional_moments):
    Vector_Operator(),
    include_double_dimensional_moments_(include_double_dimensional_moments),
    dimensional_moments_(spatial_discretization->dimensional_moments()),
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
    
    local_number_of_dimensional_moments_ = (include_double_dimensional_moments
                                            ? dimensional_moments_->number_of_double_dimensional_moments()
                                            : dimensional_moments_->number_of_dimensional_moments());
    row_size_ = psi_size;
    column_size_ = phi_size * local_number_of_dimensional_moments_;
    
    check_class_invariants();
}

void SUPG_Moment_To_Discrete::
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

    #pragma omp parallel for schedule(dynamic, 10)
    for (int i = 0; i < number_of_points; ++i)
    {
        double tau = spatial_discretization_->weight(i)->options()->tau;
        for (int g = 0; g < number_of_groups; ++g)
        {
            for (int n = 0; n < number_of_nodes; ++n)
            {
                for (int o = 0; o < number_of_ordinates; ++o)
                {
                    // Get dimensional coefficients
                    vector<double> const direction = angular_discretization_->direction(o);
                    vector<double> const coefficients
                        = (include_double_dimensional_moments_
                           ? dimensional_moments_->double_coefficients(tau,
                                                                       direction)
                           : dimensional_moments_->coefficients(tau,
                                                                direction));
                    
                    // Sum over all moments
                    double dir_sum = 0;
                    for (int m = 0; m < number_of_moments; ++m)
                    {
                        // Get dimensional summation for this moment and ordinate
                        double dim_sum = 0;
                        for (int d = 0; d < local_number_of_dimensional_moments_; ++d)
                        {
                            int k = n + number_of_nodes * (d + local_number_of_dimensional_moments_ * (g + number_of_groups * (m + number_of_moments * i)));
                            dim_sum += coefficients[d] * y[k];
                        }
                            
                        // Apply spherical harmonic function
                        int l = scattering_indices[m];
                        double p = angular_discretization_->moment(m, o);
                        dir_sum += (2 * static_cast<double>(l) + 1) / angular_normalization * p * dim_sum;
                    }

                    // Set the new discrete flux
                    int k = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                    x[k] = dir_sum;
                }
            }
        }
    }
}

void SUPG_Moment_To_Discrete::
check_class_invariants() const
{
    Assert(spatial_discretization_);
    Assert(angular_discretization_);
    Assert(energy_discretization_);
}
