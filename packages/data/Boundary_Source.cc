#include "Boundary_Source.hh"

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "XML_Node.hh"

using namespace std;

Boundary_Source::
Boundary_Source(int index,
                shared_ptr<Angular_Discretization> angular_discretization,
                shared_ptr<Energy_Discretization> energy_discretization,
                vector<double> const &boundary_source,
                vector<double> const &alpha):
    index_(index),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    boundary_source_(boundary_source),
    alpha_(alpha)
{
    check_class_invariants();
}

void Boundary_Source::
check_class_invariants() const
{
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_groups = energy_discretization_->number_of_groups();

    Assert(angular_discretization_);
    Assert(energy_discretization_);
    Assert(alpha_.size() == number_of_groups);
    Assert(boundary_source_.size() == number_of_ordinates * number_of_groups);
}

bool Boundary_Source::
has_reflection() const
{
    int number_of_groups = energy_discretization_->number_of_groups();
    
    for (int i = 0; i < number_of_groups; ++i)
    {
        if (alpha_[i] != 0)
        {
            return true;
        }
    }
    
    return false;
}

void Boundary_Source::
output(XML_Node output_node) const
{
    output_node.set_attribute(index_, "index");
    output_node.set_child_vector(alpha_, "alpha", "boundary_cell");
    output_node.set_child_vector(boundary_source_, "boundary_source", "group-ordinate");
}


