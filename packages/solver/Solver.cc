#include "Solver.hh"

#include <iomanip>
#include <iostream>
#include <string>

using namespace std;

Solver::
Solver(int solver_print):
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
