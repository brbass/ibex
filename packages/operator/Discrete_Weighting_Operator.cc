#include "Discrete_Weighting_Operator.hh"

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weight_Function.hh"

using std::shared_ptr;
using std::vector;

Discrete_Weighting_Operator::
Discrete_Weighting_Operator(shared_ptr<Weak_Spatial_Discretization> spatial,
                            shared_ptr<Angular_Discretization> angular,
                            shared_ptr<Energy_Discretization> energy):
    Square_Vector_Operator(spatial->number_of_points()
                           * spatial->number_of_nodes()
                           * angular->number_of_ordinates()
                           * energy->number_of_groups()),
    spatial_(spatial),
    angular_(angular),
    energy_(energy)
{
}

void Discrete_Weighting_Operator::
apply(vector<double> &x) const
{
    // Get size data
    int number_of_points = spatial_->number_of_points();
    int number_of_nodes = spatial_->number_of_nodes();
    int number_of_groups = energy_->number_of_groups();
    int number_of_ordinates = angular_->number_of_ordinates();
    int dimension = spatial_->dimension();
    int number_of_dimensional_moments = spatial_->number_of_dimensional_moments();
    
    vector<double> result(number_of_points * number_of_nodes * number_of_groups * number_of_ordinates, 0);
    
    for (int i = 0; i < number_of_points; ++i)
    {
        // Get weight function and data
        shared_ptr<Weight_Function> weight = spatial_->weight(i);
        int number_of_basis_functions = weight->number_of_basis_functions();
        vector<int> basis_indices = weight->basis_function_indices();
        Weight_Function::Material_Options options = weight->material_options();
        double const tau = options.tau;
        bool const normalized = options.normalized;
        bool const include_supg = options.include_supg;
        vector<double> const iv_w = weight->iv_w();
        vector<double> const iv_dw = weight->iv_dw();
        vector<double> const iv_b_w = weight->iv_b_w();
        vector<double> const iv_b_dw = weight->iv_b_dw();
        
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            // Get normalization constant
            vector<double> const direction = angular_->direction(o);
            double norm = 1;
            if (!normalized)
            {
                norm = iv_w[0];
                if (include_supg)
                {
                    for (int d = 0; d < dimension; ++d)
                    {
                        norm += tau * direction[d] * iv_dw[d];
                    }
                }
            }
            
            for (int j = 0; j < number_of_basis_functions; ++j)
            {
                // Get summation constant
                double mult = iv_b_w[j];
                if (include_supg)
                {
                    for (int d = 0; d < dimension; ++d)
                    {
                        int k = d + dimension * j;
                        mult += tau * direction[d] * iv_b_dw[k];
                    }
                }
                mult /= norm;
                
                for (int g = 0; g < number_of_groups; ++g)
                {
                    for (int n = 0; n < number_of_nodes; ++n)
                    {
                        int k_bas = basis_indices[j];
                        int k_x = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * k_bas));
                        int k_res = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                        
                        result[k_res] += mult * x[k_x];
                    }
                }
            }
        }
    }
    
    // Put result into "x"
    result.swap(x);
}

