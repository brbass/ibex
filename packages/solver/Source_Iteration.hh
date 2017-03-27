#ifndef Source_Iteration_hh
#define Source_Iteration_hh

#include "Solver.hh"

class Source_Iteration : public Solver
{
public:

    struct Options
    {
        int max_iterations = 1000;
        int solver_print = 0;
        double tolerance = 1e-8;
    };
    
    struct Result
    {
        int source_iterations;
        int total_iterations;
        std::vector<double> phi;
    };
    
    Source_Iteration(Options options);

    virtual void solve() override;

    virtual void get_flux(std::vector<double> &x) const override;

    virtual void output(XML_Node output_node) const override;
    
    virtual void check_class_invariants() const override;
    
private:
    
    Options options_;
    Result result_;
};

#endif
