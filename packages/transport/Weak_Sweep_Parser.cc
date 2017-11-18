#include "Weak_Sweep_Parser.hh"

#include "Angular_Discretization.hh"
#include "Conversion.hh"
#include "Energy_Discretization.hh"
#include "Transport_Discretization.hh"
#include "Weak_RBF_Sweep.hh"
#include "Weak_Spatial_Discretization.hh"
#include "XML_Node.hh"

using namespace std;

Weak_Sweep_Parser::
Weak_Sweep_Parser(shared_ptr<Weak_Spatial_Discretization> spatial,
                  shared_ptr<Angular_Discretization> angular,
                  shared_ptr<Energy_Discretization> energy,
                  shared_ptr<Transport_Discretization> transport):
    spatial_(spatial),
    angular_(angular),
    energy_(energy),
    transport_(transport)
{
}

shared_ptr<Weak_RBF_Sweep> Weak_Sweep_Parser::
get_weak_rbf_sweep(XML_Node input_node) const
{
    Weak_RBF_Sweep::Options options;
    options.quit_if_diverged = input_node.get_attribute<bool>("quit_if_diverged",
                                                              options.quit_if_diverged);
    options.max_iterations = input_node.get_attribute<int>("max_iterations",
                                                           options.max_iterations);
    options.max_restarts = input_node.get_attribute<int>("max_restarts",
                                                         options.max_restarts);
    options.kspace = input_node.get_attribute<int>("kspace",
                                                   options.kspace);
    options.level_of_fill = input_node.get_attribute<double>("level_of_fill",
                                                             options.level_of_fill);
    options.tolerance = input_node.get_attribute<double>("tolerance",
                                                         options.tolerance);
    options.drop_tolerance = input_node.get_attribute<double>("drop_tolerance",
                                                              options.drop_tolerance);
    
    string solver = input_node.get_attribute<string>("solver",
                                                     "amesos");
    options.solver = options.solver_conversion()->convert(solver);
    
    return make_shared<Weak_RBF_Sweep>(options,
                                       spatial_,
                                       angular_,
                                       energy_,
                                       transport_);
}
