#include "Energy_Discretization.hh"

#include <vector>

#include "Check.hh"
#include "XML_Node.hh"

using namespace std;

Energy_Discretization::
Energy_Discretization(int number_of_groups):
    number_of_groups_(number_of_groups),
    energy_bounds_(number_of_groups + 1, 0)
{
    check_class_invariants();
}

Energy_Discretization::
Energy_Discretization(int number_of_groups,
                      vector<double> const &energy_bounds):
    number_of_groups_(number_of_groups),
    energy_bounds_(energy_bounds)
{
    check_class_invariants();
}

void Energy_Discretization::
check_class_invariants() const
{
    Assert(energy_bounds_.size() == number_of_groups_ + 1);
}

void Energy_Discretization::
output(XML_Node output_node) const
{
    output_node.set_child_value(number_of_groups_, "number_of_groups");
    output_node.set_child_vector(energy_bounds_, "energy_bounds", "group-bound");
}
