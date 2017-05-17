#include "Transport_Discretization.hh"

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Spatial_Discretization.hh"
#include "XML_Node.hh"

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
    has_reflection_ = spatial->has_reflection();
    
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

void Transport_Discretization::
output(XML_Node output_node) const
{
    output_node.set_attribute(has_reflection_, "has_reflection");
    output_node.set_child_value(phi_size_, "phi_size");
    output_node.set_child_value(psi_size_, "psi_size");
    output_node.set_child_value(number_of_augments_, "number_of_augments");
}
