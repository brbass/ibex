#include "Dimensional_Moment_Copy.hh"

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weight_Function.hh"

using std::shared_ptr;
using std::vector;

Dimensional_Moment_Copy::
Dimensional_Moment_Copy(shared_ptr<Weak_Spatial_Discretization> spatial,
                        shared_ptr<Angular_Discretization> angular,
                        shared_ptr<Energy_Discretization> energy):
    Vector_Operator(spatial->number_of_points()
                    * spatial->number_of_nodes()
                    * energy->number_of_groups()
                    * angular->number_of_moments(),
                    spatial->number_of_points()
                    * spatial->number_of_nodes()
                    * spatial->number_of_dimensional_moments()
                    * energy->number_of_groups()
                    * angular->number_of_moments()),
                    
    spatial_(spatial),
    angular_(angular),
    energy_(energy)
{
}

void Dimensional_Moment_Copy::
check_class_invariants() const
{
    Assert(spatial_);
    Assert(angular_);
    Assert(energy_);
    Assert(spatial_->number_of_dimensional_moments() == spatial_->dimension() + 1);
}

void Dimensional_Moment_Copy::
apply(vector<double> &x) const
{
    // Get data
    int number_of_points = spatial_->number_of_points();
    int number_of_nodes = spatial_->number_of_nodes();
    int number_of_dimensional_moments = spatial_->number_of_dimensional_moments();
    int number_of_groups = energy_->number_of_groups();
    int number_of_moments = angular_->number_of_moments();

    // Resize data
    int size = number_of_points * number_of_nodes * number_of_groups * number_of_moments;
    vector<double> result(size * number_of_dimensional_moments);

    // Put data into result
    for (int i = 0; i < size; ++i)
    {
        for (int d = 0; d < number_of_dimensional_moments; ++d)
        {
            result[d + size * i] = x[i];
        }
    }

    // Swap result with x
    x.swap(result);
}
