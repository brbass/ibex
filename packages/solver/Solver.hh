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
    
    struct Result
    {
        // General
        int total_iterations = -1;
        std::vector<double> coefficients;
        std::vector<std::vector<double> > phi;

        // Source iteration
        int source_iterations = -1;

        // Eigenvalue
        double k_eigenvalue = -1;
    };
    
    // Constructor
    Solver(int solver_print,
           Solver::Type type);
    
    // Solve problem
    virtual void solve() = 0;
    
    // Ouput data to XML file
    virtual void output(XML_Node output_node) const = 0;
    virtual void output_result(XML_Node output_node,
                               std::shared_ptr<Result> result) const;
    
    // Check class values
    virtual void check_class_invariants() const = 0;

    virtual std::shared_ptr<Result> result() const = 0;
    
protected:
    
    // Print functions
    virtual void print_name(std::string solution_type) const final;
    virtual void print_iteration(int iteration) const final;
    virtual void print_convergence() const final;
    virtual void print_value(double value) const final;
    virtual void print_error(double error) const final;
    virtual void print_failure() const final;
    virtual void print_eigenvalue(double eigenvalue) const final;

    // Solver print code
    int solver_print_;
    Solver::Type type_;
};

#endif
