#include "Solver_Parser.hh"

#include "Identity_Operator.hh"
#include "Krylov_Eigenvalue.hh"
#include "Krylov_Steady_State.hh"
#include "Linf_Convergence.hh"
#include "Moment_Value_Operator.hh"
#include "Solver_Factory.hh"
#include "Source_Iteration.hh"
#include "Vector_Operator.hh"
#include "Vector_Operator_Functions.hh"
#include "Weak_Spatial_Discretization.hh"
#include "XML_Node.hh"

using namespace std;

Solver_Parser::
Solver_Parser(shared_ptr<Weak_Spatial_Discretization> spatial,
              shared_ptr<Angular_Discretization> angular,
              shared_ptr<Energy_Discretization> energy,
              shared_ptr<Transport_Discretization> transport):
    spatial_(spatial),
    angular_(angular),
    energy_(energy),
    transport_(transport),
    factory_(make_shared<Solver_Factory>(spatial,
                                         angular,
                                         energy,
                                         transport))
{
}

shared_ptr<Source_Iteration> Solver_Parser::
get_source_iteration(XML_Node input_node,
                     shared_ptr<Sweep_Operator> Linv) const
{
    // Get combined operators
    shared_ptr<Vector_Operator> source_operator;
    shared_ptr<Vector_Operator> flux_operator;
    factory_->get_source_operators(Linv,
                                   source_operator,
                                   flux_operator);

    // Get value operator
    vector<shared_ptr<Vector_Operator> > value_operators
        = {make_shared<Moment_Value_Operator>(spatial_,
                                              angular_,
                                              energy_,
                                              false)};

    // Get convergence
    shared_ptr<Convergence_Measure> convergence
        = make_shared<Linf_Convergence>();

    
    // Get options
    Source_Iteration::Options iteration_options;
    iteration_options.max_source_iterations = input_node.get_attribute<int>("max_source_iterations", 5000);
    iteration_options.max_iterations = input_node.get_attribute<int>("max_iterations", 5000);
    iteration_options.solver_print = input_node.get_attribute<int>("solver_print", 0);
    iteration_options.tolerance = input_node.get_attribute<double>("tolerance", 1e-10);
    
    // Create solver
    return make_shared<Source_Iteration>(iteration_options,
                                         spatial_,
                                         angular_,
                                         energy_,
                                         transport_,
                                         convergence,
                                         source_operator,
                                         flux_operator,
                                         value_operators); 
}

shared_ptr<Krylov_Steady_State> Solver_Parser::
get_krylov_steady_state(XML_Node input_node,
                        shared_ptr<Sweep_Operator> Linv) const
{
    // Get combined operators
    shared_ptr<Vector_Operator> source_operator;
    shared_ptr<Vector_Operator> flux_operator;
    factory_->get_source_operators(Linv,
                                   source_operator,
                                   flux_operator);
    shared_ptr<Identity_Operator> identity
        = make_shared<Identity_Operator>(flux_operator->column_size());
    flux_operator = identity - flux_operator;
    
    // Get value operator
    vector<shared_ptr<Vector_Operator> > value_operators
        = {make_shared<Moment_Value_Operator>(spatial_,
                                              angular_,
                                              energy_,
                                              false)};

    // Get convergence
    shared_ptr<Convergence_Measure> convergence
        = make_shared<Linf_Convergence>();

    // Get options
    Krylov_Steady_State::Options iteration_options;
    iteration_options.max_source_iterations = input_node.get_attribute<int>("max_source_iterations", 5000);
    iteration_options.max_iterations = input_node.get_attribute<int>("max_iterations", 1000);
    iteration_options.kspace = input_node.get_attribute<int>("kspace", 10);
    iteration_options.solver_print = input_node.get_attribute<int>("solver_print", 0);
    iteration_options.tolerance = input_node.get_attribute<double>("tolerance", 1e-10);
    
    // Create solver
    return make_shared<Krylov_Steady_State>(iteration_options,
                                            spatial_,
                                            angular_,
                                            energy_,
                                            transport_,
                                            convergence,
                                            source_operator,
                                            flux_operator,
                                            value_operators); 
}

std::shared_ptr<Krylov_Eigenvalue> Solver_Parser::
get_krylov_eigenvalue(XML_Node input_node,
                      shared_ptr<Sweep_Operator> Linv) const
{
    // Get combined operators
    shared_ptr<Vector_Operator> fission_operator;
    shared_ptr<Vector_Operator> flux_operator;
    factory_->get_eigenvalue_operators(Linv,
                                       fission_operator,
                                       flux_operator);
    shared_ptr<Identity_Operator> identity
        = make_shared<Identity_Operator>(flux_operator->column_size());
    flux_operator = identity - flux_operator;
    
    // Get value operator
    vector<shared_ptr<Vector_Operator> > value_operators
        = {make_shared<Moment_Value_Operator>(spatial_,
                                              angular_,
                                              energy_,
                                              false)};

    // Get source iteration
    Krylov_Eigenvalue::Options iteration_options;
    iteration_options.max_inverse_iterations = input_node.get_attribute<int>("max_inverse_iterations", 5000);
    iteration_options.max_iterations = input_node.get_attribute<int>("max_iterations", 1000);
    iteration_options.kspace = input_node.get_attribute<int>("kspace", 10);
    iteration_options.solver_print = input_node.get_attribute<int>("solver_print", 0);
    iteration_options.tolerance = input_node.get_attribute<double>("tolerance", 1e-10);
    
    return make_shared<Krylov_Eigenvalue>(iteration_options,
                                          spatial_,
                                          angular_,
                                          energy_,
                                          transport_,
                                          fission_operator,
                                          flux_operator,
                                          value_operators); 
}
