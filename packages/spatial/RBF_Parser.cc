#include "RBF_Parser.hh"

#include <string>

#include "Check.hh"
#include "Compact_Gaussian_RBF.hh"
#include "Gaussian_RBF.hh"
#include "Inverse_Multiquadric_RBF.hh"
#include "Multiquadric_RBF.hh"
#include "RBF.hh"
#include "Truncated_Gaussian_RBF.hh"
#include "Wendland1_RBF.hh"
#include "Wendland3_RBF.hh"
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
    if (rbf_type == "gaussian")
    {
        return make_shared<Gaussian_RBF>();
    }
    else if (rbf_type == "multiquadric")
    {
        return make_shared<Multiquadric_RBF>();
    }
    else if (rbf_type == "inverse_multiquadric")
    {
        return make_shared<Inverse_Multiquadric_RBF>();
    }
    else if (rbf_type == "compact_gaussian")
    {
        double radius = input_node.get_attribute<double>("radius",
                                                         5.0);
        return make_shared<Compact_Gaussian_RBF>(radius);
    }
    else if (rbf_type == "truncated_gaussian")
    {
        double radius = input_node.get_attribute<double>("radius",
                                                         5.0);
        return make_shared<Truncated_Gaussian_RBF>(radius);
    }
    else if (rbf_type == "wendland10")
    {
        return make_shared<Wendland1_RBF>(0);
    }
    else if (rbf_type == "wendland11")
    {
        return make_shared<Wendland1_RBF>(1);
    }
    else if (rbf_type == "wendland12")
    {
        return make_shared<Wendland1_RBF>(2);
    }
    else if (rbf_type == "wendland13")
    {
        return make_shared<Wendland1_RBF>(3);
    }
    else if (rbf_type == "wendland30")
    {
        return make_shared<Wendland3_RBF>(0);
    }
    else if (rbf_type == "wendland31")
    {
        return make_shared<Wendland3_RBF>(1);
    }
    else if (rbf_type == "wendland32")
    {
        return make_shared<Wendland3_RBF>(2);
    }
    else if (rbf_type == "wendland33")
    {
        return make_shared<Wendland3_RBF>(3);
    }
    else
    {
        AssertMsg(false, "basis_type \"" + rbf_type + "\" not found");
        return shared_ptr<RBF>();
    }
}
