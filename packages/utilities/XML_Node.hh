#ifndef XML_Node_hh
#define XML_Node_hh

#include <memory>
#include <string>

#include "pugixml.hh"

/*
  Interface for pugi::xml_node class
*/
class XML_Node
{
public:

    // Find a child node
    XML_Node get_child(string name);

    // Append a child node
    XML_Node append_child(string name);
    
protected:

    // Create an XML_Node (see XML_Document for public creation)
    XML_Node(pugi::xml_node node);
    
private:
    
    pugi::xml_node xml_node_;
};

#endif
