#include "Arbitrary_Moment_Value_Operator.hh"

#if defined(ENABLE_OPENMP)
    #include <omp.h>
#endif

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weight_Function.hh"

using std::shared_ptr;
using std::vector;

Arbitrary_Moment_Value_Operator::
Arbitrary_Moment_Value_Operator(shared_ptr<Weak_Spatial_Discretization> spatial,
                                shared_ptr<Angular_Discretization> angular,
                                shared_ptr<Energy_Discretization> energy,
                                vector<vector<double> > const &evaluation_points):
    Vector_Operator(),
    spatial_(spatial),
    angular_(angular),
    energy_(energy),
    number_of_evaluation_points_(evaluation_points.size()),
    evaluation_points_(evaluation_points)
{
    column_size_ = (spatial->number_of_points()
                    * spatial->number_of_nodes()
                    * angular->number_of_moments()
                    * energy->number_of_groups());
    row_size_ = (number_of_evaluation_points_
                    * spatial->number_of_nodes()
                    * angular->number_of_moments()
                    * energy->number_of_groups());
    check_class_invariants();
}

void Arbitrary_Moment_Value_Operator::
apply(vector<double> &x) const
{
    // Get size data
    int number_of_points = spatial_->number_of_points();
    int number_of_nodes = spatial_->number_of_nodes();
    int number_of_groups = energy_->number_of_groups();
    int number_of_moments = angular_->number_of_moments();
    int number_of_values = number_of_nodes * number_of_groups * number_of_moments;
    
    // Initialize result
    vector<double> result(number_of_evaluation_points_ * number_of_nodes * number_of_groups * number_of_moments, 0);
    #pragma omp parallel for schedule(dynamic, 10)
    for (int i = 0; i < number_of_evaluation_points_; ++i)
    {
        // Get expansion values
        vector<double> values
            = spatial_->expansion_values(number_of_values,
                                         evaluation_points_[i],
                                         x);
        
        // Assign expansion values to result
        for (int m = 0; m < number_of_moments; ++m)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    int k_val = n + number_of_nodes * (g + number_of_groups * m);
                    int k_res = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    
                    result[k_res] = values[k_val];
                }
            }
        }
    }
    
    // Put result into "x"
    x.swap(result);
}

void Arbitrary_Moment_Value_Operator::
check_class_invariants() const
{
}
