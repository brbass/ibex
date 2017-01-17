#include "Gauss_Legendre_Quadrature.hh"

#include "Check.hh"
#include "Quadrature_Rule.hh"
#include "XML_Node.hh"

using std::vector;

Gauss_Legendre_Quadrature::
Gauss_Legendre_Quadrature(int dimension,
                          int number_of_moments,
                          int number_of_ordinates):
    Angular_Discretization(dimension,
                           number_of_moments,
                           number_of_ordinates)
{
    Quadrature_Rule::gauss_legendre(number_of_ordinates, ordinates_, weights_);
    
    directions_.resize(number_of_ordinates_);
    
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        vector<double> direction = {ordinates_[o]};

        directions_[o] = direction;
    }
    
    check_class_invariants();
}

void Gauss_Legendre_Quadrature::
check_class_invariants() const
{
    Assert(dimension_ == 1);
    Assert(ordinates_.size() == number_of_ordinates_);
    Assert(weights_.size() == number_of_ordinates_);
    Assert(directions_.size() == number_of_ordinates_);
}

void Gauss_Legendre_Quadrature::
output(XML_Node output_node) const
{
    output_node.set_attribute("gauss_legendre_quadrature", "quadrature_type");
    output_node.set_child_value(dimension_, "dimension");
    output_node.set_child_value(number_of_moments_, "number_of_moments");
    output_node.set_child_value(number_of_ordinates_, "number_of_ordinates");
    output_node.set_child_vector(ordinates_, "ordinates", "ordinate");
    output_node.set_child_vector(weights_, "weights", "ordinate");
}
