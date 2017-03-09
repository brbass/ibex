#include "Boundary_Source.hh"

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "XML_Node.hh"

using namespace std;

Boundary_Source::
Boundary_Source(int index,
                Dependencies dependencies,
                shared_ptr<Angular_Discretization> angular_discretization,
                shared_ptr<Energy_Discretization> energy_discretization,
                vector<double> const &boundary_source,
                vector<double> const &alpha):
    index_(index),
    dependencies_(dependencies),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    alpha_(alpha)
{
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    // int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_groups = energy_discretization_->number_of_groups();
    size_ = number_of_ordinates * number_of_groups;
    
    switch (dependencies_.angular)
    {
    case Dependencies::Angular::ISOTROPIC:
        Assert(boundary_source.size() == number_of_groups); 
        boundary_source_.resize(size_);
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                int index = g + number_of_groups * o;
                boundary_source_[index] = boundary_source[g];
            }
        }
        break;
    case Dependencies::Angular::MOMENTS:
        AssertMsg(false, "not yet implemented");
        break;
    case Dependencies::Angular::ORDINATES:
        Assert(boundary_source.size() == size_); 
        boundary_source_ = boundary_source;
        break;
    }
    
    check_class_invariants();
}

void Boundary_Source::
check_class_invariants() const
{
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_groups = energy_discretization_->number_of_groups();

    Assert(angular_discretization_);
    Assert(energy_discretization_);
    Assert(alpha_.size() == number_of_groups);
    Assert(boundary_source_.size() == size_);
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


