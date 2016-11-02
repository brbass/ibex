#ifndef Solver_hh
#define Solver_hh

#include <memory>
#include <string>
#include <vector>

#include "pugixml.hh"

using std::shared_ptr;
using std::string;
using std::vector;

/*
  Pure virtual class for a solver of the transport equation
*/
class Solver
{
public:
    
    // Constructor
    Solver(int solver_print);
    
    // Solve fixed source problem
    virtual void solve_steady_state(vector<double> &x) = 0;
    
    // Solve k-eigenvalue problem
    virtual void solve_k_eigenvalue(double &k_eigenvalue, vector<double> &x) = 0;
    
    // Solve time-dependent problem
    virtual void solve_time_dependent(vector<double> &x) = 0;
    
    // Ouput data to XML file
    virtual void output(pugi::xml_node &output_node) const = 0;

protected:

    virtual void print_name(string solution_type) const;
    virtual void print_iteration(int iteration) const;
    virtual void print_convergence() const;
    virtual void print_value(double value) const;
    virtual void print_error(double error) const;
    virtual void print_failure() const;
    virtual void print_eigenvalue(double eigenvalue) const;
    
    int solver_print_;
};

#endif
