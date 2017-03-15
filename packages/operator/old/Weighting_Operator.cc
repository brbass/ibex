#include "Weighting_Operator.hh"

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weight_Function.hh"

using std::shared_ptr;
using std::vector;

Weighting_Operator::
Weighting_Operator(shared_ptr<Weak_Spatial_Discretization> spatial,
                   shared_ptr<Angular_Discretization> angular,
                   shared_ptr<Energy_Discretization> energy):
    Vector_Operator(spatial->number_of_points()
                    * spatial->number_of_nodes()
                    * angular->number_of_moments()
                    * energy->number_of_groups(),
                    spatial->number_of_points()
                    * spatial->number_of_nodes()
                    * spatial->number_of_dimensional_moments()
                    * angular->number_of_moments()
                    * energy->number_of_groups()),
    spatial_(spatial),
    angular_(angular),
    energy_(energy)
{
}
    
void Weighting_Operator::
apply(vector<double> &x) const
{
    // Get size data
    int number_of_points = spatial_->number_of_points();
    int number_of_nodes = spatial_->number_of_nodes();
    int number_of_groups = energy_->number_of_groups();
    int number_of_moments = angular_->number_of_moments();
    int number_of_dimensional_moments = spatial_->number_of_dimensional_moments();
    int dimension = spatial_->dimension();

    // Initialize result
    vector<double> result(number_of_points * number_of_nodes * number_of_groups * number_of_moments * number_of_dimensional_moments);
    
    for (int i = 0; i < number_of_points; ++i)
    {
        // Get weight function and integral information
        shared_ptr<Weight_Function> weight = spatial_->weight(i);
        int number_of_basis_functions = weight->number_of_basis_functions();
        vector<int> basis_indices = weight->basis_function_indices();
        vector<double> const iv_b_w = weight->iv_b_w();
        vector<double> const iv_b_dw = weight->iv_b_dw();

        // Each moment, group and node is done independently
        for (int m = 0; m < number_of_moments; ++m)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    // Perform sum over standard integral
                    {
                        double sum = 0;
                        for (int j = 0; j < number_of_basis_functions; ++j)
                        {
                            int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * basis_indices[j]));
                            sum += iv_b_w[j] * x[k_phi];
                        }
                        int k_res = n + number_of_nodes * (0 + number_of_dimensional_moments * (g + number_of_groups * (m + number_of_moments * i)));
                        result[k_res] = sum;
                    }

                    // Perform sum over gradient integrals
                    for (int d = 1; d < number_of_dimensional_moments; ++d)
                    {
                        double sum = 0;
                        for (int j = 0; j < number_of_basis_functions; ++j)
                        {
                            int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * basis_indices[j]));
                            int k_int = (d - 1) + dimension * j;
                            sum += iv_b_dw[k_int] * x[k_phi];
                        }
                        int k_res = n + number_of_nodes * (d + number_of_dimensional_moments * (g + number_of_groups * (m + number_of_moments * i)));
                        result[k_res] = sum;
                    }
                }
            }
        }
    }
}
