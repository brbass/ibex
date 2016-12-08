#include "XML_Node.hh"

using namespace std;

XML_Node::
XML_Node(pugi::xml_node xml_node):
    xml_node_(xml_node)
{
}

XML_Node XML_Node::
get_child(string name,
          bool check)
{
    pugi::xml_node child_node = xml_node_.child(name.c_str());
    
    if (child_node.empty() && check)
    {
        string parent_name = static_cast<std::string>(xml_node_.name());
        AssertMsg(false, "child node (" + name + ") in node (" + parent_name + ") not found");
    }
    
    return XML_Node(child_node);
}

XML_Node XML_Node::
get_sibling(string name,
            bool check)
{
    pugi::xml_node sibling_node = xml_node_.next_sibling(name.c_str());

    if (sibling_node.empty() && check)
    {
        string self_name = static_cast<std::string>(xml_node_.name());
        AssertMsg(false, "sibling node (" + name + ") in node (" + self_name + ") not found");
    }
    
    return sibling_node;
}

XML_Node XML_Node::
append_child(string name)
{
    return XML_Node(xml_node_.append_child(name.c_str()));
}
