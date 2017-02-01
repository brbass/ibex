#include "RBF_Parser.hh"

#include <string>

#include "RBF.hh"
#include "XML_Node.hh"

using namespace std;

RBF_Parser::
RBF_Parser()
{
}

shared_ptr<RBF> RBF_Parser::
parse_from_xml(XML_Node input_node)
{
    string rbf_type = input_node.get_attribute<string>("function");
    
    if (rbf_type == "compact_gaussian")
    {
        double radius = input_node.get_child_value<double>("radius",
                                                              5.0);
        return make_shared<Compact_Gaussian_RBF>(radius);
    }
    else if (rbf_type == "truncated_gaussian")
    {
        double radius = input_node.get_child_value<double>("radius",
                                                              5.0);
        return make_shared<Truncated_Gaussian_RBF>(radius);
    }
    else if (rbf_type == "wendland30")
    {
        return make_shared<Wendland_RBF>(0);
    }
    else if (rbf_type == "wendland31")
    {
        return make_shared<Wendland_RBF>(1);
    }
    else if (rbf_type == "wendland32")
    {
        return make_shared<Wendland_RBF>(2);
    }
    else if (rbf_type == "wendland33")
    {
        return make_shared<Wendland_RBF>(3);
    }
    else
    {
        AssertMsg(false, "basis_type \"" + rbf_type + "\" not found");
        return shared_ptr<RBF>();
    }
}
