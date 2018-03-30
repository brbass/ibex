#include "Meshless_Sweep_Parser.hh"

#include "Angular_Discretization.hh"
#include "Conversion.hh"
#include "Energy_Discretization.hh"
#include "Strong_Meshless_Sweep.hh"
#include "Transport_Discretization.hh"
#include "Weak_Meshless_Sweep.hh"
#include "Weak_Spatial_Discretization.hh"
#include "XML_Node.hh"

using namespace std;

Meshless_Sweep_Parser::
Meshless_Sweep_Parser(shared_ptr<Weak_Spatial_Discretization> spatial,
                      shared_ptr<Angular_Discretization> angular,
                      shared_ptr<Energy_Discretization> energy,
                      shared_ptr<Transport_Discretization> transport):
    spatial_(spatial),
    angular_(angular),
    energy_(energy),
    transport_(transport)
{
}

shared_ptr<Meshless_Sweep> Meshless_Sweep_Parser::
get_meshless_sweep(XML_Node input_node) const
{
    switch (spatial_->options()->discretization)
    {
    case Weak_Spatial_Discretization_Options::Discretization::WEAK:
        return get_weak_sweep(input_node);
    case Weak_Spatial_Discretization_Options::Discretization::STRONG:
        return get_strong_sweep(input_node);
    }
    
}

shared_ptr<Meshless_Sweep> Meshless_Sweep_Parser::
get_weak_sweep(XML_Node input_node) const
{
    Meshless_Sweep::Options options;
    options.quit_if_diverged = input_node.get_attribute<bool>("quit_if_diverged",
                                                              options.quit_if_diverged);
    options.use_preconditioner = input_node.get_attribute<bool>("use_preconditioner",
                                                                options.use_preconditioner);
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
    options.weighted_preconditioner = input_node.get_attribute<bool>("weighted_preconditioner",
                                                                     options.weighted_preconditioner);
    
    string solver = input_node.get_attribute<string>("solver",
                                                     "amesos");
    options.solver = options.solver_conversion()->convert(solver);
    
    return make_shared<Weak_Meshless_Sweep>(options,
                                       spatial_,
                                       angular_,
                                       energy_,
                                       transport_);
}

shared_ptr<Meshless_Sweep> Meshless_Sweep_Parser::
get_strong_sweep(XML_Node input_node) const
{
    Meshless_Sweep::Options options;
    options.quit_if_diverged = input_node.get_attribute<bool>("quit_if_diverged",
                                                              options.quit_if_diverged);
    options.use_preconditioner = input_node.get_attribute<bool>("use_preconditioner",
                                                                options.use_preconditioner);
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
    
    return make_shared<Strong_Meshless_Sweep>(options,
                                              spatial_,
                                              angular_,
                                              energy_,
                                              transport_);
}
