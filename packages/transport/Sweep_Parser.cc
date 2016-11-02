#include "Sweep_Parser.hh"

#include "Check.hh"
#include "Discrete_To_Moment.hh"
#include "Fission.hh"
#include "Moment_To_Discrete.hh"
#include "Preconditioner.hh"
#include "RBF_Collocation_Sweep.hh"
#include "Sweep_Operator.hh"
#include "Scattering.hh"
#include "Transport_Discretization.hh"

using namespace std;

Sweep_Parser::
Sweep_Parser(pugi::xml_node &input_file,
             shared_ptr<Spatial_Discretization> spatial,
             shared_ptr<Angular_Discretization> angular,
             shared_ptr<Energy_Discretization> energy,
             shared_ptr<Transport_Discretization> transport):
    Parser(input_file),
    spatial_(spatial),
    angular_(angular),
    energy_(energy),
    transport_(transport)
{
    sweeper_ = parse_sweeper();
}    

shared_ptr<Sweep_Operator> Sweep_Parser::
parse_sweeper()
{
    pugi::xml_node sweeper = input_file_.child("sweep_operator");
    
    string sweeper_type = XML_Functions::child_value<string>(sweeper, "type");
    
    if (sweeper_type == "rbf_collocation")
    {
        return parse_rbf_collocation(sweeper);
    }
    else
    {
        AssertMsg(false, "sweeper type " + sweeper_type + " not found");
        
        return shared_ptr<Sweep_Operator>();
    }
}

shared_ptr<RBF_Collocation_Sweep> Sweep_Parser::
parse_rbf_collocation(pugi::xml_node &input_node)
{
    shared_ptr<RBF_Discretization> rbf_discretization
        = dynamic_pointer_cast<RBF_Discretization>(spatial_);

    Assert(rbf_discretization);

    string solution_variable_str = XML_Functions::child_value<string>(input_node, "solution_variable");
    
    RBF_Collocation_Sweep::Solution_Variable solution_variable;
    
    if (solution_variable_str == "coefficient")
    {
        solution_variable = RBF_Collocation_Sweep::Solution_Variable::COEFFICIENT;
    }
    else if (solution_variable_str == "psi")
    {
        solution_variable = RBF_Collocation_Sweep::Solution_Variable::PSI;
    }
    else
    {
        AssertMsg(false, "solution variable " + solution_variable_str + " not found");
    }

    string matrix_str = XML_Functions::child_value<string>(input_node, "matrix_solver");
    RBF_Collocation_Sweep::Matrix_Solver matrix_solver;

    if (matrix_str == "amesos")
    {
        matrix_solver = RBF_Collocation_Sweep::Matrix_Solver::AMESOS;
    }
    else if (matrix_str == "aztec")
    {
        matrix_solver = RBF_Collocation_Sweep::Matrix_Solver::AZTEC;
    }
    else
    {
        AssertMsg(false, "matrix solver " + matrix_str + " not found");
    }
    
    string condition_str = XML_Functions::child_value<string>(input_node, "condition_calculation");
    RBF_Collocation_Sweep::Condition_Calculation condition_calculation;
    
    if (condition_str == "cheap")
    {
        condition_calculation = RBF_Collocation_Sweep::Condition_Calculation::CHEAP;
    }
    else if (condition_str == "expensive")
    {
        condition_calculation = RBF_Collocation_Sweep::Condition_Calculation::EXPENSIVE;
    }
    else if (condition_str == "none")
    {
        condition_calculation = RBF_Collocation_Sweep::Condition_Calculation::NONE;
    }
    else if (condition_str == "aztec")
    {
        condition_calculation = RBF_Collocation_Sweep::Condition_Calculation::AZTEC;
    }
    else
    {
        AssertMsg(false, "condition calculation " + condition_str + " not found");
    }
    
    return make_shared<RBF_Collocation_Sweep>(solution_variable,
                                              matrix_solver,
                                              condition_calculation,
                                              rbf_discretization,
                                              angular_,
                                              energy_,
                                              transport_);
}
