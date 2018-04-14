#include "Discrete_Value_Operator.hh"

#if defined(ENABLE_OPENMP)
    #include <omp.h>
#endif

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weight_Function.hh"

using std::shared_ptr;
using std::vector;

Discrete_Value_Operator::
Discrete_Value_Operator(shared_ptr<Weak_Spatial_Discretization> spatial,
                        shared_ptr<Angular_Discretization> angular,
                        shared_ptr<Energy_Discretization> energy,
                        bool weighted):
    Square_Vector_Operator(),
    spatial_(spatial),
    angular_(angular),
    energy_(energy),
    weighted_(weighted)
{
    size_ = (spatial->number_of_points()
             * spatial->number_of_nodes()
             * angular->number_of_ordinates()
             * energy->number_of_groups());
    
    check_class_invariants();
}

void Discrete_Value_Operator::
apply(vector<double> &x) const
{
    // Get size data
    int number_of_points = spatial_->number_of_points();
    int number_of_nodes = spatial_->number_of_nodes();
    int number_of_groups = energy_->number_of_groups();
    int number_of_ordinates = angular_->number_of_ordinates();

    vector<double> y(x);
    x.assign(number_of_points * number_of_nodes * number_of_groups * number_of_ordinates, 0);

    #pragma omp parallel for schedule(dynamic, 10)
    for (int i = 0; i < number_of_points; ++i)
    {
        // Get weight function and data
        shared_ptr<Weight_Function> weight = spatial_->weight(i);
        int number_of_basis_functions = weight->number_of_basis_functions();
        Weight_Function::Integrals const integrals = weight->integrals();
        Weight_Function::Values const values = weight->values();
        vector<int> basis_indices = weight->basis_function_indices();
        double const iv_w = (weighted_
                             ? integrals.iv_w[0]
                             : 1);
        vector<double> const &iv_b_w = (weighted_
                                        ? integrals.iv_b_w
                                        : values.v_b);
        
        for (int j = 0; j < number_of_basis_functions; ++j)
        {
            // Get summation constant
            double const mult = iv_b_w[j] / iv_w;
            int k_bas = basis_indices[j];
            
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                for (int g = 0; g < number_of_groups; ++g)
                {
                    for (int n = 0; n < number_of_nodes; ++n)
                    {
                        int k_y = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * k_bas));
                        int k_x = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                        
                        x[k_x] += mult * x[k_y];
                    }
                }
            }
        }
    }
}

void Discrete_Value_Operator::
check_class_invariants() const
{
    Assert(spatial_);
    Assert(angular_);
    Assert(energy_);
}
