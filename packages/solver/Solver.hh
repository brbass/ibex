#ifndef Solver_hh
#define Solver_hh

#include <memory>
#include <string>
#include <vector>

class XML_Node;

/*
  Pure virtual class for a solver of the transport equation
*/
class Solver
{
public:
    
    enum class Type
    {
        STEADY_STATE,
        K_EIGENVALUE,
        ALPHA_EIGENVALUE,
        TIME_DEPENDENT
    };
    
    // Constructor
    Solver(int solver_print,
           Solver::Type type);
    
    // Solve problem
    virtual void solve() = 0;
    
    // Get solution
    virtual void get_eigenvalue(double &eigenvalue) const;
    virtual void get_flux(std::vector<double> &x) const = 0;
    
    // Ouput data to XML file
    virtual void output(XML_Node output_node) const = 0;

    // Check class values
    virtual void check_class_invariants() const = 0;
    
protected:
    
    // Print functions
    virtual void print_name(std::string solution_type) const;
    virtual void print_iteration(int iteration) const;
    virtual void print_convergence() const;
    virtual void print_value(double value) const;
    virtual void print_error(double error) const;
    virtual void print_failure() const;
    virtual void print_eigenvalue(double eigenvalue) const;

    // Solver print code
    int solver_print_;
    Solver::Type type_;
};

#endif
