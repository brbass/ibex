#include "RBF_Factory.hh"

#include <string>

#include "Compact_Gaussian_RBF.hh"
#include "Gaussian_RBF.hh"
#include "Inverse_Multiquadric_RBF.hh"
#include "Multiquadric_RBF.hh"
#include "RBF.hh"
#include "Truncated_Gaussian_RBF.hh"
#include "Wendland_RBF.hh"

using std::make_shared;
using std::shared_ptr;
using std::string;

RBF_Factory::
RBF_Factory()
{
}

shared_ptr<RBF> RBF_Factory::
get_rbf(string rbf_type)
{
    if (rbf_type == "gaussian")
    {
        return make_shared<Gaussian_RBF>();
    }
    else if (rbf_type == "multiquadric_rbf")
    {
        return make_shared<Multiquadric_RBF>();
    }
    else if (rbf_type == "inverse_multiquadric_rbf")
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
