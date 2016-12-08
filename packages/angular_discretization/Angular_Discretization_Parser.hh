#ifndef Angular_Discretization_Parser_hh
#define Angular_Discretization_Parser_hh

#include <memory>

class Angular_Discretization;
class XML_Node;

/* 
   Create an object of type Angular_Discretization from xml input file
*/
class Angular_Discretization_Parser
{
public:
    
    // Creator
    Angular_Discretization_Parser();
    
    // Parse from XML node
    std::shared_ptr<Angular_Discretization> parse_from_xml(XML_Node input_node);
};

#endif
