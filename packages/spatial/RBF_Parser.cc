#include "RBF_Parser.hh"

#include <string>

#include "RBF.hh"
#include "RBF_Factory.hh"
#include "XML_Node.hh"

using std::make_shared;
using std::shared_ptr;
using std::string;

RBF_Parser::
RBF_Parser()
{
}

shared_ptr<RBF> RBF_Parser::
parse_from_xml(XML_Node input_node)
{
    string rbf_type = input_node.get_attribute<string>("function");
    RBF_Factory factory;
    
    return factory.get_rbf(rbf_type);
}
