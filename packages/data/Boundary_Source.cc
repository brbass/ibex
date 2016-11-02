#include "Boundary_Source.hh"

#include <cmath>
#include <string>

#include "Check.hh"
#include "XML_Functions.hh"

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
output(pugi::xml_node &output_node) const
{
    pugi::xml_node source = output_node.append_child("boundary_source");
    
    XML_Functions::append_child(source, alpha_, "alpha", "boundary_cell");
    XML_Functions::append_child(source, boundary_source_, "boundary_source", "group-ordinate");
}


