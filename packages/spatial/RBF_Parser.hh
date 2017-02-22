#ifndef RBF_Parser_hh
#define RBF_Parser_hh

#include <memory>

class RBF;
class XML_Node;

class RBF_Parser
{
public:

    RBF_Parser();
    
    std::shared_ptr<RBF> parse_from_xml(XML_Node input_node);
};

#endif
