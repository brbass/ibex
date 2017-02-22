#ifndef Energy_Discretization_Parser_hh
#define Energy_Discretization_Parser_hh

#include <memory>

class Energy_Discretization;
class XML_Node;

/*
  Create an Energy_Discretization object from XML input file
*/
class Energy_Discretization_Parser
{
public:

    // Constructor
    Energy_Discretization_Parser();

    // Parse from XML node
    std::shared_ptr<Energy_Discretization> parse_from_xml(XML_Node input_node);
};

#endif
