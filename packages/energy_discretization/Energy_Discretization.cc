#include "Energy_Discretization.hh"

#include <vector>

#include "Check.hh"
#include "XML_Functions.hh"

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
output(pugi::xml_node &output_node) const
{
    pugi::xml_node energy = output_node.append_child("energy_discretization");
    
    XML_Functions::append_child(energy, number_of_groups_, "number_of_groups");
    XML_Functions::append_child(energy, energy_bounds_, "energy_bounds", "group_bound");
}
