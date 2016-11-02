#include "Energy_Discretization_Parser.hh"

#include "XML_Functions.hh"

using namespace std;

Energy_Discretization_Parser::
Energy_Discretization_Parser(pugi::xml_node &input_file):
    Parser(input_file)
{
    pugi::xml_node energy = input_file.child("energy_discretization");
    
    int number_of_groups = XML_Functions::child_value<int>(energy, "number_of_groups");
    
    energy_ = make_shared<Energy_Discretization>(number_of_groups);
}
