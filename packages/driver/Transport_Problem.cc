#include "Transport_Problem.hh"

#include "Solver.hh"
#include "XML_Functions.hh"

Transport_Problem::
Transport_Problem(Problem_Type problem_type,
                  shared_ptr<Solver> solver):
    problem_type_(problem_type),
    solver_(solver)
{
}

void Transport_Problem::
solve()
{
    switch(problem_type_)
    {
    case Problem_Type::STEADY_STATE:
        solver_->solve_steady_state(phi_);
        break;
    case Problem_Type::K_EIGENVALUE:
        solver_->solve_k_eigenvalue(k_eigenvalue_, phi_);
        break;
    case Problem_Type::TIME_DEPENDENT:
        solver_->solve_time_dependent(phi_);
        break;
    }
}

void Transport_Problem::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node solution = output_node.append_child("solution");

    if (problem_type_ == Problem_Type::K_EIGENVALUE)
    {
        XML_Functions::append_child(solution, k_eigenvalue_, "k_eigenvalue");
    }
    XML_Functions::append_child(solution, phi_, "phi", "node-group-moment-cell");
    // XML_Functions::append_child(solution, psi_, "psi");
}
