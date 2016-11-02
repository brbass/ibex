#include "Solver_Parser.hh"

#include "Check.hh"
#include "Discrete_To_Moment.hh"
#include "Fission.hh"
#include "Krylov_Iteration.hh"
#include "Moment_To_Discrete.hh"
#include "Null_Solver.hh"
#include "Power_Iteration.hh"
#include "Preconditioner.hh"
#include "Sweep_Operator.hh"
#include "Sweep_Parser.hh"
#include "Scattering.hh"
#include "Source_Iteration.hh"
#include "Transport_Discretization.hh"

using namespace std;

Solver_Parser::
Solver_Parser(pugi::xml_node &input_file,
              shared_ptr<Spatial_Discretization> spatial,
              shared_ptr<Angular_Discretization> angular,
              shared_ptr<Energy_Discretization> energy,
              shared_ptr<Transport_Discretization> transport,
              shared_ptr<Sweep_Operator> sweeper):
    Parser(input_file),
    spatial_(spatial),
    angular_(angular),
    energy_(energy),
    transport_(transport),
    sweeper_(sweeper)
{
    pugi::xml_node solver_node = input_file.child("solution_method");
    
    string solver_type = XML_Functions::child_value<string>(solver_node, "type");
    
    if (solver_type == "source_iteration")
    {
        solver_  = parse_source_iteration();
    }
    else if (solver_type == "krylov_iteration")
    {
        solver_ = parse_krylov_iteration();
    }
    else if (solver_type == "power_iteration")
    {
        solver_ = parse_power_iteration();
    }
    else if (solver_type == "null_solver")
    {
        solver_ = make_shared<Null_Solver>(1,
                                           spatial,
                                           angular,
                                           energy);
    }
    else
    {
        AssertMsg(false, "solver type " + solver_type + " not found");
    }
}    

shared_ptr<Source_Iteration> Solver_Parser::
parse_source_iteration()
{
    pugi::xml_node solver_node = input_file_.child("solution_method");
    
    int max_iterations = XML_Functions::child_value<int>(solver_node, "max_iterations");
    int solver_print = XML_Functions::child_value<int>(solver_node, "solver_print");
    double tolerance = XML_Functions::child_value<double>(solver_node, "tolerance");
    shared_ptr<Discrete_To_Moment> discrete_to_moment = parse_discrete_to_moment();
    shared_ptr<Moment_To_Discrete> moment_to_discrete = parse_moment_to_discrete();
    shared_ptr<Scattering> scattering = parse_scattering();
    shared_ptr<Fission> fission = parse_fission();
    shared_ptr<Preconditioner> preconditioner = parse_preconditioner(solver_node,
                                                                     scattering,
                                                                     fission);
    
    return make_shared<Source_Iteration>(max_iterations,
                                         solver_print,
                                         tolerance,
                                         spatial_,
                                         angular_,
                                         energy_,
                                         transport_,
                                         sweeper_,
                                         discrete_to_moment,
                                         moment_to_discrete,
                                         scattering,
                                         fission,
                                         preconditioner);
}

shared_ptr<Krylov_Iteration> Solver_Parser::
parse_krylov_iteration()
{
    pugi::xml_node solver_node = input_file_.child("solution_method");
    
    int max_iterations = XML_Functions::child_value<int>(solver_node, "max_iterations");
    int kspace = XML_Functions::child_value<int>(solver_node, "kspace");
    int solver_print = XML_Functions::child_value<int>(solver_node, "solver_print");
    double tolerance = XML_Functions::child_value<double>(solver_node, "tolerance");
    shared_ptr<Discrete_To_Moment> discrete_to_moment = parse_discrete_to_moment();
    shared_ptr<Moment_To_Discrete> moment_to_discrete = parse_moment_to_discrete();
    shared_ptr<Scattering> scattering = parse_scattering();
    shared_ptr<Fission> fission = parse_fission();
    shared_ptr<Preconditioner> preconditioner = parse_preconditioner(solver_node,
                                                                     scattering,
                                                                     fission);
    
    return make_shared<Krylov_Iteration>(max_iterations,
                                         kspace,
                                         solver_print,
                                         tolerance,
                                         spatial_,
                                         angular_,
                                         energy_,
                                         transport_,
                                         sweeper_,
                                         discrete_to_moment,
                                         moment_to_discrete,
                                         scattering,
                                         fission,
                                         preconditioner);
}

shared_ptr<Power_Iteration> Solver_Parser::
parse_power_iteration()
{
    pugi::xml_node solver_node = input_file_.child("solution_method");
    
    int max_iterations = XML_Functions::child_value<int>(solver_node, "max_iterations");
    int kspace = XML_Functions::child_value<int>(solver_node, "kspace");
    int solver_print = XML_Functions::child_value<int>(solver_node, "solver_print");
    double tolerance = XML_Functions::child_value<double>(solver_node, "tolerance");
    shared_ptr<Discrete_To_Moment> discrete_to_moment = parse_discrete_to_moment();
    shared_ptr<Moment_To_Discrete> moment_to_discrete = parse_moment_to_discrete();
    shared_ptr<Scattering> scattering = parse_scattering();
    shared_ptr<Fission> fission = parse_fission();
    shared_ptr<Preconditioner> preconditioner = parse_preconditioner(solver_node,
                                                                      scattering,
                                                                      fission);
    
    return make_shared<Power_Iteration>(max_iterations,
                                        kspace,
                                        solver_print,
                                        tolerance,
                                        spatial_,
                                        angular_,
                                        energy_,
                                        transport_,
                                        sweeper_,
                                        discrete_to_moment,
                                        moment_to_discrete,
                                        scattering,
                                        fission,
                                        preconditioner);
}


shared_ptr<Discrete_To_Moment> Solver_Parser::
parse_discrete_to_moment()
{
    return make_shared<Discrete_To_Moment>(spatial_,
                                           angular_,
                                           energy_);
}

shared_ptr<Moment_To_Discrete> Solver_Parser::
parse_moment_to_discrete()
{
    return make_shared<Moment_To_Discrete>(spatial_,
                                           angular_,
                                           energy_);
}

shared_ptr<Scattering> Solver_Parser::
parse_scattering()
{
    return make_shared<Scattering>(spatial_,
                                   angular_,
                                   energy_);
}

shared_ptr<Fission> Solver_Parser::
parse_fission()
{
    return make_shared<Fission>(spatial_,
                                angular_,
                                energy_);
}

shared_ptr<Preconditioner> Solver_Parser::
parse_preconditioner(pugi::xml_node solver_node,
                     shared_ptr<Scattering_Operator> scattering,
                     shared_ptr<Scattering_Operator> fission)
{
    return shared_ptr<Preconditioner>();
}
