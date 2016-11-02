#include "Angular_Discretization_Parser.hh"

#include "Gauss_Legendre_Quadrature.hh"
#include "LDFE_Quadrature.hh"
#include "XML_Functions.hh"

using namespace std;

Angular_Discretization_Parser::
Angular_Discretization_Parser(pugi::xml_node &input_file):
    Parser(input_file)
{
    pugi::xml_node angular = input_file.child("angular_discretization");
    
    int dimension = XML_Functions::child_value<int>(angular, "dimension");
    int number_of_moments = XML_Functions::child_value<int>(angular, "number_of_moments");

    if (dimension == 1)
    {
        int number_of_ordinates = XML_Functions::child_value<int>(angular, "number_of_ordinates");

        angular_ = make_shared<Gauss_Legendre_Quadrature>(dimension,
                                                          number_of_moments,
                                                          number_of_ordinates);
    }
    else
    {
        int rule = XML_Functions::child_value<int>(angular, "rule");
        
        angular_ = make_shared<LDFE_Quadrature>(dimension,
                                                number_of_moments,
                                                rule);
    }
}
