#ifndef Source_Iteration_hh
#define Source_Iteration_hh

#include "Solver.hh"

class Angular_Discretization;
class Convergence_Measure;
class Energy_Discretization;
class Spatial_Discretization;
class Transport_Discretization;
class Vector_Operator;

class Source_Iteration : public Solver
{
public:

    struct Options
    {
        int max_source_iterations = 1000;
        int max_iterations = 1000;
        int solver_print = 0;
        double tolerance = 1e-8;
    };
    
    struct Result
    {
        int source_iterations = -1;
        int total_iterations = -1;
        std::vector<double> coefficients;
        std::vector<std::vector<double> > phi;
    };
    
    Source_Iteration(Options options,
                     std::shared_ptr<Spatial_Discretization> spatial_discretization,
                     std::shared_ptr<Angular_Discretization> angular_discretization,
                     std::shared_ptr<Energy_Discretization> energy_discretization,
                     std::shared_ptr<Transport_Discretization> transport_discretization,
                     std::shared_ptr<Convergence_Measure> convergence,
                     std::shared_ptr<Vector_Operator> source_operator,
                     std::shared_ptr<Vector_Operator> flux_operator,
                     std::vector<std::shared_ptr<Vector_Operator> > value_operators);
    
    virtual void solve() override;
    std::shared_ptr<Result> result() const
    {
        return result_;
    }
    virtual void output(XML_Node output_node) const override;
    
    virtual void check_class_invariants() const override;
    
private:

    // Input data
    Options options_;
    std::shared_ptr<Spatial_Discretization> spatial_discretization_;
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;
    std::shared_ptr<Transport_Discretization> transport_discretization_;
    std::shared_ptr<Convergence_Measure> convergence_;
    std::shared_ptr<Vector_Operator> source_operator_;
    std::shared_ptr<Vector_Operator> flux_operator_;
    std::vector<std::shared_ptr<Vector_Operator> > value_operators_;
    
    
    // Output data
    std::shared_ptr<Result> result_;
};

#endif
