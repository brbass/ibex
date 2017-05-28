#include "Discrete_Normalization_Operator.hh"

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weight_Function.hh"

using std::shared_ptr;
using std::vector;

Discrete_Normalization_Operator::
Discrete_Normalization_Operator(shared_ptr<Weak_Spatial_Discretization> spatial,
                                shared_ptr<Angular_Discretization> angular,
                                shared_ptr<Energy_Discretization> energy,
                                Options options):
    Weighting_Operator(spatial,
                       angular,
                       energy,
                       options)
{
    size_ = (spatial->number_of_points()
             * spatial->number_of_nodes()
             * angular->number_of_ordinates()
             * energy->number_of_groups());
    
    check_class_invariants();
}

void Discrete_Normalization_Operator::
apply(vector<double> &x) const
{
    // Get size data
    int number_of_points = spatial_->number_of_points();
    int number_of_nodes = spatial_->number_of_nodes();
    int number_of_groups = energy_->number_of_groups();
    int number_of_ordinates = angular_->number_of_ordinates();
    int dimension = spatial_->dimension();
    int number_of_dimensional_moments = spatial_->number_of_dimensional_moments();
    shared_ptr<Weak_Spatial_Discretization_Options> const weak_options
        = spatial_->weak_options();
    
    bool include_normalization;
    switch (options_.normalization)
    {
    case Options::Normalization::AUTO:
        include_normalization = !weak_options->normalized;
        break;
    case Options::Normalization::TRUE:
        include_normalization = true;
        break;
    case Options::Normalization::FALSE:
        include_normalization = false;
        break;
    }
        
    bool include_supg;
    switch (options_.include_supg)
    {
    case Options::Include_SUPG::AUTO:
        include_supg = weak_options->include_supg;
        break;
    case Options::Include_SUPG::TRUE:
        include_supg = true;
        break;
    case Options::Include_SUPG::FALSE:
        include_supg = false;
        break;
    }
        
    if (include_normalization)
    {
        for (int i = 0; i < number_of_points; ++i)
        {
            // Get weight function and data
            shared_ptr<Weight_Function> weight = spatial_->weight(i);
            int number_of_basis_functions = weight->number_of_basis_functions();
            vector<int> basis_indices = weight->basis_function_indices();
            Weight_Function::Integrals const integrals = weight->integrals();
            shared_ptr<Weight_Function_Options> const weight_options = weight->options();
            double const tau = weight_options->tau;
            vector<double> const &iv_w = integrals.iv_w;
            vector<double> const &iv_dw = integrals.iv_dw;
        
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                // Get normalization constant
                vector<double> const direction = angular_->direction(o);
                double norm = 1;
                norm = iv_w[0];
                if (include_supg)
                {
                    for (int d = 0; d < dimension; ++d)
                    {
                        norm += tau * direction[d] * iv_dw[d];
                    }
                }
                double invnorm = 1. / norm;
                // Apply normalization constant
                for (int g = 0; g < number_of_groups; ++g)
                {
                    for (int n = 0; n < number_of_nodes; ++n)
                    {
                        int k = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                    
                        x[k] *= invnorm;
                    }
                }
            }
        }
    }
}

void Discrete_Normalization_Operator::
check_class_invariants() const
{
    Assert(spatial_);
    Assert(angular_);
    Assert(energy_);
}
