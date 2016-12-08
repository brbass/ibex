#include "Energy_Discretization_Parser.hh"

#include "Energy_Discretization.hh"
#include "XML_Node.hh"

using namespace std;

Energy_Discretization_Parser::
Energy_Discretization_Parser()
{
}

shared_ptr<Energy_Discretization> Energy_Discretization_Parser::
parse_from_xml(XML_Node input_node)
{
    int number_of_groups = input_node.get_child_value<int>("number_of_groups");
    
    return make_shared<Energy_Discretization>(number_of_groups);
}
