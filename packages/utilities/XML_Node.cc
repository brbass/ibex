#include "XML_Node.hh"

using namespace std;

XML_Node::
XML_Node(pugi::xml_node xml_node):
    xml_node_(xml_node)
{
}

XML_Node XML_Node::
get_child(string name)
{
    return XML_Node(xml_node_.child(name.c_str()));
}

XML_Node XML_Node::
append_child(string name)
{
    return XML_Node(xml_node_.append_child(name.c_str()));
}
