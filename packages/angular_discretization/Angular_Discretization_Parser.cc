#include "Angular_Discretization_Parser.hh"

#include "Angular_Discretization.hh"
#include "Gauss_Legendre_Quadrature.hh"
#include "LDFE_Quadrature.hh"
#include "XML_Node.hh"

using namespace std;

Angular_Discretization_Parser::
Angular_Discretization_Parser()
{
}

shared_ptr<Angular_Discretization> Angular_Discretization_Parser::
parse_from_xml(XML_Node input_node)
{
    int dimension = input_node.get_child_value<int>("dimension");
    int number_of_moments = input_node.get_child_value<int>("number_of_moments");
    
    if (dimension == 1)
    {
        int number_of_ordinates = input_node.get_child_value<int>("number_of_ordinates");
        
        return make_shared<Gauss_Legendre_Quadrature>(dimension,
                                                      number_of_moments,
                                                      number_of_ordinates);
    }
    else
    {
        int rule = input_node.get_child_value<int>("rule");
        
        return make_shared<LDFE_Quadrature>(dimension,
                                            number_of_moments,
                                            rule);
    }
}
