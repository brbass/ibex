#include "Moment_Weighting_Operator.hh"

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weight_Function.hh"

using std::shared_ptr;
using std::vector;

Moment_Weighting_Operator::
Moment_Weighting_Operator(shared_ptr<Weak_Spatial_Discretization> spatial,
                          shared_ptr<Angular_Discretization> angular,
                          shared_ptr<Energy_Discretization> energy,
                          Options options):
    Weighting_Operator(spatial,
                       angular,
                       energy,
                       options)
{
    int phi_size = (spatial->number_of_points()
                    * spatial->number_of_nodes()
                    * angular->number_of_moments()
                    * energy->number_of_groups());
    int number_of_dimensional_moments
        = spatial->number_of_dimensional_moments();
    
    switch (options.include_supg)
    {
    case Options::Include_SUPG::AUTO:
        // fallthrough to true
    case Options::Include_SUPG::TRUE:
        local_number_of_dimensional_moments_ = number_of_dimensional_moments;
        break;
    case Options::Include_SUPG::FALSE:
        local_number_of_dimensional_moments_ = 1;
        break;
    }
    
    column_size_ = phi_size;
    row_size_ = phi_size * local_number_of_dimensional_moments_;
    
    check_class_invariants();
}

void Moment_Weighting_Operator::
apply(vector<double> &x) const
{
    // Get size data
    int number_of_points = spatial_->number_of_points();
    int number_of_nodes = spatial_->number_of_nodes();
    int number_of_groups = energy_->number_of_groups();
    int number_of_moments = angular_->number_of_moments();
    int dimension = spatial_->dimension();
    
    vector<double> result(number_of_points * number_of_nodes * number_of_groups * number_of_moments * local_number_of_dimensional_moments_, 0);
    
    for (int i = 0; i < number_of_points; ++i)
    {
        // Get weight function and data
        shared_ptr<Weight_Function> weight = spatial_->weight(i);
        int number_of_basis_functions = weight->number_of_basis_functions();
        Weight_Function::Integrals const integrals = weight->integrals();
        Weight_Function::Values const values = weight->values();
        vector<int> basis_indices = weight->basis_function_indices();
        Weight_Function::Options weight_options = weight->options();
        bool include_normalization;
        switch (options_.normalization)
        {
        case Options::Normalization::AUTO:
            if (weight_options.normalized)
            {
                include_normalization = false;
            }
            else
            {
                if (weight_options.include_supg)
                {
                    include_normalization = false;
                }
                else
                {
                    include_normalization = true;
                }
            }
            break;
        case Options::Normalization::TRUE:
            include_normalization = true;
            break;
        case Options::Normalization::FALSE:
            include_normalization = false;
            break;
        }

        vector<double> const &iv_w = integrals.iv_w;
        vector<double> const &iv_b_w = integrals.iv_b_w;
        vector<double> const &iv_b_dw = integrals.iv_b_dw;
        
        // Get normalization constant
        double const norm = include_normalization ? iv_w[0] : 1;
        
        for (int j = 0; j < number_of_basis_functions; ++j)
        {
            // Get summation constants and other basis data
            vector<double> mult(local_number_of_dimensional_moments_);
            mult[0] = iv_b_w[j] / norm;
            for (int d = 1; d < local_number_of_dimensional_moments_; ++d)
            {
                int k_i = d - 1 + dimension * j;
                mult[d] = iv_b_dw[k_i] / norm;
            }
            int k_bas = basis_indices[j];

            // Apply weighting to each moment, group, node and (if SUPG) dimensional moment
            for (int m = 0; m < number_of_moments; ++m)
            {
                for (int g = 0; g < number_of_groups; ++g)
                {
                    for (int n = 0; n < number_of_nodes; ++n)
                    {
                        int k_x = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * k_bas));
                        for (int d = 0; d < local_number_of_dimensional_moments_; ++d)
                        {
                            int k_res = n + number_of_nodes * (d + local_number_of_dimensional_moments_ * (g + number_of_groups * (m + number_of_moments * i)));
                            
                            result[k_res] += mult[d] * x[k_x];
                        }
                    }
                }
            }
        }
    }
    
    // Put result into "x"
    x.swap(result);
}

void Moment_Weighting_Operator::
check_class_invariants() const
{
    Assert(spatial_);
    Assert(angular_);
    Assert(energy_);
    if (options_.include_supg == Options::Include_SUPG::TRUE)
    {
        Assert(options_.normalization == Options::Normalization::FALSE);
    }
}
