#include "XML_Node.hh"

using namespace std;

XML_Node::
XML_Node(shared_ptr<pugi::xml_node> xml_node,
         string name):
    xml_node_(xml_node),
    name_(name)
{
}

XML_Node XML_Node::
get_child(string name,
          bool check)
{
    shared_ptr<pugi::xml_node> child_node = make_shared<pugi::xml_node>(xml_node_->child(name.c_str()));
    
    if (child_node->empty() && check)
    {
        AssertMsg(false, "child node (" + name + ") in node (" + name_ + ") not found");
    }

    string child_name = name_ + "/" + name;
    int name_size = child_name.size() - 32;
    if (name_size > 0)
    {
        child_name = "..." + child_name.substr(name_size);
    }
    
    return XML_Node(child_node,
                    child_name);
}

XML_Node XML_Node::
get_sibling(string name,
            bool check)
{
    shared_ptr<pugi::xml_node> sibling_node = make_shared<pugi::xml_node>(xml_node_->next_sibling(name.c_str()));
    
    if (sibling_node->empty() && check)
    {
        AssertMsg(false, "sibling node (" + name + ") in node (" + name_ + ") not found");
    }
    
    string sibling_name = name_ + "/" + name;
    int name_size = sibling_name.size() - 32;
    if (name_size > 0)
    {
        sibling_name = "..." + sibling_name.substr(name_size);
    }
    
    return XML_Node(sibling_node,
                    sibling_name);
}

XML_Node XML_Node::
append_child(string name)
{
    string child_name = name_ + "/" + name;
    int name_size = child_name.size() - 32;
    if (name_size > 0)
    {
        child_name = "..." + child_name.substr(name_size);
    }
    
    return XML_Node(make_shared<pugi::xml_node>(xml_node_->append_child(name.c_str())),
                    child_name);
}

void XML_Node::
prepend_all(XML_Node copy_node)
{
    shared_ptr<pugi::xml_node> copy_pugi_node
        = copy_node.xml_node();
    
    for (pugi::xml_node node = copy_pugi_node->first_child();
         node;
         node = node.next_sibling())
    {
        xml_node_->prepend_copy(node);
    }
}

void XML_Node::
append_all(XML_Node copy_node)
{
    shared_ptr<pugi::xml_node> copy_pugi_node
        = copy_node.xml_node();
    
    for (pugi::xml_node node = copy_pugi_node->first_child();
         node;
         node = node.next_sibling())
    {
        xml_node_->append_copy(node);
    }
}

void XML_Node::
prepend_node(XML_Node copy_node)
{
    xml_node_->prepend_copy(*copy_node.xml_node())
}

void XML_Node::
append_node(XML_Node copy_node)
{
    xml_node_->append_copy(*copy_node.xml_node())
}
