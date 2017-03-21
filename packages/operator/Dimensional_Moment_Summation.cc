#include "Dimensional_Moment_Summation.hh"

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weight_Function.hh"

using std::shared_ptr;
using std::vector;

Dimensional_Moment_Summation::
Dimensional_Moment_Summation(shared_ptr<Weak_Spatial_Discretization> spatial,
                             shared_ptr<Angular_Discretization> angular,
                             shared_ptr<Energy_Discretization> energy):
    Vector_Operator(),                    
    spatial_(spatial),
    angular_(angular),
    energy_(energy)
{
    int phi_size =(spatial->number_of_points()
                   * spatial->number_of_nodes()
                   * energy->number_of_groups()
                   * angular->number_of_moments());
    row_size_ = phi_size * spatial->number_of_dimensional_moments();
    column_size_ = phi_size;
}

void Dimensional_Moment_Summation::
check_class_invariants() const
{
    Assert(spatial_);
    Assert(angular_);
    Assert(energy_);
    Assert(spatial_->number_of_dimensional_moments() == spatial_->dimension() + 1);
}

void Dimensional_Moment_Summation::
apply(vector<double> &x) const
{
    // Get data
    int number_of_points = spatial_->number_of_points();
    int number_of_nodes = spatial_->number_of_nodes();
    int number_of_dimensional_moments = spatial_->number_of_dimensional_moments();
    int number_of_groups = energy_->number_of_groups();
    int number_of_ordinates = angular_->number_of_ordinates();

    // Initialize vectors for coefficients and result
    vector<double> coefficients(number_of_dimensional_moments);
    coefficients[0] = 1;
    vector<double> result(number_of_points * number_of_nodes * number_of_groups * number_of_ordinates);

    // Sum the dimensional coefficients
    for (int i = 0; i < number_of_points; ++i)
    {
        double tau = spatial_->weight(i)->options().tau;
        
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            vector<double> const direction = angular_->direction(o);
            for (int d = 1; d < number_of_dimensional_moments; ++d)
            {
                coefficients[d] = direction[d - 1] * tau;
            }
            
            for (int g = 0; g < number_of_groups; ++g)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    double sum = 0;
                    for (int d = 0; d < number_of_dimensional_moments; ++d)
                    {
                        int k_psi_d = n + number_of_nodes * (d + number_of_dimensional_moments * (g + number_of_groups * (o + number_of_ordinates * i)));
                        sum += x[k_psi_d] * coefficients[d];
                    }
                    
                    int k_psi = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                    result[k_psi] = sum;
                }
            }
        }
    }
    
    x.swap(result);
}
