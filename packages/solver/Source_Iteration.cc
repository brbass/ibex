#include "Source_Iteration.hh"

#include "XML_Node.hh"

using namespace std;

Source_Iteration::
Source_Iteration(Options options):
    Solver(options.solver_print,
           Solver::Type::STEADY_STATE),
    options_(options)
{
}

void Source_Iteration::
solve()
{
    
}

void Source_Iteration::
get_flux(vector<double> &x) const
{
    
}

void Source_Iteration::
output(XML_Node output_node) const
{
    
}

void Source_Iteration::
check_class_invariants() const
{
    
}
