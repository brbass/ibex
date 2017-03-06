#include "Discrete_To_Moment.hh"

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "Spatial_Discretization.hh"

using std::shared_ptr;
using std::vector;

Discrete_To_Moment::
Discrete_To_Moment(shared_ptr<Spatial_Discretization> spatial_discretization,
                   shared_ptr<Angular_Discretization> angular_discretization,
                   shared_ptr<Energy_Discretization> energy_discretization,
                   bool include_dimensional_moments):
    Vector_Operator(spatial_discretization->number_of_points()
                    * spatial_discretization->number_of_nodes()
                    * (include_dimensional_moments
                       ? spatial_discretization->number_of_dimensional_moments()
                       : 1)
                    * energy_discretization->number_of_groups()
                    * angular_discretization->number_of_moments(), // row size
                    spatial_discretization->number_of_points()
                    * spatial_discretization->number_of_nodes()
                    * (include_dimensional_moments
                       ? spatial_discretization->number_of_dimensional_moments()
                       : 1)
                    * energy_discretization->number_of_groups()
                    * angular_discretization->number_of_ordinates()), // column size
    include_dimensional_moments_(include_dimensional_moments),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization)
{
    check_class_invariants();
}

void Discrete_To_Moment::
apply(vector<double> &x) const
{
    vector<double> y(x);
    
    x.resize(row_size());
    
    int number_of_points = spatial_discretization_->number_of_points();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int local_number_of_dimensional_moments = (include_dimensional_moments_
                                               ? spatial_discretization_->number_of_dimensional_moments()
                                               : 1);
    vector<double> const weights = angular_discretization_->weights();
    vector<double> const ordinates = angular_discretization_->ordinates();
    
    for (int i = 0; i < number_of_points; ++i)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            for (int d = 0; d < local_number_of_dimensional_moments; ++d)
            {
            for (int n = 0; n < number_of_nodes; ++n)
            {
                for (int m = 0; m < number_of_moments; ++m)
                {
                    double sum = 0;
                    
                    for (int o = 0; o < number_of_ordinates; ++o)
                    {
                        int k = n + number_of_nodes * (d + local_number_of_dimensional_moments * (g + number_of_groups * (o + number_of_ordinates * i)));
                        
                        double p = angular_discretization_->moment(m, o);

                        sum += weights[o] * p * y[k];
                    }
                    
                    int k = n + number_of_nodes * (d + local_number_of_dimensional_moments * (g + number_of_groups * (m + number_of_moments * i)));
                    
                    x[k] = sum;
                }
            }
            }
        }
    }
}

void Discrete_To_Moment::
check_class_invariants() const
{
    Assert(spatial_discretization_);
    Assert(angular_discretization_);
    Assert(energy_discretization_);
}
