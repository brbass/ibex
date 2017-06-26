#include "Solver.hh"

#include <iomanip>
#include <iostream>
#include <string>

#include "Check.hh"
#include "XML_Node.hh"

using namespace std;

Solver::
Solver(int solver_print,
       Solver::Type type):
    solver_print_(solver_print)
{
}

void Solver::
print_name(string solution_type) const
{
    if (solver_print_)
    {
        cout << endl;
        cout << "\t\t*******************************************************";
        cout << endl;
        cout << "\t\t***** " << solution_type;
        cout << endl;
        cout << "\t\t*******************************************************";
        cout << endl;
    }        
}

void Solver::
print_iteration(int iteration) const
{
    if (solver_print_)
    {
        cout << "\t\titer:\t";
        cout << iteration;
        cout << "\t";
    }
}

void Solver::
print_convergence() const
{
    if (solver_print_)
    {
        cout << endl;
        cout << "\t\tConverged";
        cout << endl;
    }
}

void Solver::
print_failure() const
{
    if (solver_print_)
    {
        cout << endl;
        cout << "\t\tFailed to converge";
        cout << endl;
    }
}

void Solver::
print_value(double value) const
{
    if (solver_print_)
    {
        cout << "value:\t";
        cout << value;
        cout << "\t";
    }
}

void Solver::
print_error(double error) const
{
    if (solver_print_)
    {
        cout << "error:\t";
        cout << error;
        cout << endl;
    }
}

void Solver::
print_eigenvalue(double eigenvalue) const
{
    if (solver_print_)
    {
        cout << endl;
        cout << "\t\tk_eigenvalue:\t";
        cout << setprecision(10);
        cout << eigenvalue;
        cout << endl;
    }
}

void Solver::
output_result(XML_Node output_node,
              shared_ptr<Result> result) const
{
    output_node.set_child_value(result->total_iterations,
                                "total_iterations");
    output_node.set_child_value(result->inverse_iterations,
                                "inverse_iterations");
    output_node.set_child_vector(result->coefficients,
                                 "coefficients",
                                 "node-group-moment-point");
    XML_Node phi_node = output_node.append_child("values");
    for (vector<double> const &phi : result->phi)
    {
        phi_node.set_child_vector(phi, "phi", "node-group-moment-point");
    }

    if (result->source_iterations != -1)
    {
        output_node.set_child_value(result->source_iterations,
                                    "source_iterations");
    }
    if (result->k_eigenvalue != -1)
    {
        output_node.set_child_value(result->k_eigenvalue,
                                    "k_eigenvalue");
    }
}
