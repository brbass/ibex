#include "Transport_Discretization.hh"

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Spatial_Discretization.hh"

using std::shared_ptr;
using std::vector;

Transport_Discretization::
Transport_Discretization(shared_ptr<Spatial_Discretization> spatial,
                         shared_ptr<Angular_Discretization> angular,
                         shared_ptr<Energy_Discretization> energy):
    spatial_(spatial),
    angular_(angular),
    energy_(energy)
{
    int number_of_points = spatial->number_of_points();
    int number_of_nodes = spatial->number_of_nodes();
    int number_of_boundary_points = spatial->number_of_boundary_points();
    int number_of_groups = energy->number_of_groups();
    int number_of_ordinates = angular->number_of_ordinates();
    int number_of_moments = angular->number_of_moments();
    vector<int> const boundary_points = spatial->boundary_points();
    
    has_reflection_ = false;
    for (int i = 0; i < number_of_boundary_points; ++i)
    {
        if (spatial->point(boundary_points[i])->boundary_source()->has_reflection())
        {
            has_reflection_ = true;
            break;
        }
    }
    
    phi_size_ = number_of_points * number_of_groups * number_of_moments * number_of_nodes;
    psi_size_ = number_of_points * number_of_groups * number_of_ordinates * number_of_nodes;
    if (has_reflection_)
    {
        number_of_augments_ = number_of_boundary_points * number_of_groups * number_of_ordinates * number_of_nodes;
    }
    else
    {
        number_of_augments_ = 0;
    }
}
