#ifndef Source_Iteration_hh
#define Source_Iteration_hh

#include "Solver.hh"

class Angular_Discretization;
class Energy_Discretization;
class Spatial_Discretization;
class Transport_Discretization;
class Vector_Operator;

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
    
    Source_Iteration(Options options,
                     std::shared_ptr<Spatial_Discretization> spatial_discretization,
                     std::shared_ptr<Angular_Discretization> angular_discretization,
                     std::shared_ptr<Energy_Discretization> energy_discretization,
                     std::shared_ptr<Transport_Discretization> transport_discretization,                     std::shared_ptr<Vector_Operator> source_operator,
                     std::shared_ptr<Vector_Operator> flux_operator,
                     std::shared_ptr<Vector_Operator> value_operator);
    
    virtual void solve() override;
    
    virtual void get_flux(std::vector<double> &x) const override;

    virtual void output(XML_Node output_node) const override;
    
    virtual void check_class_invariants() const override;
    
private:

    // Input data
    Options options_;
    std::shared_ptr<Spatial_Discretization> spatial_discretization;
    std::shared_ptr<Angular_Discretization> angular_discretization;
    std::shared_ptr<Energy_Discretization> energy_discretization;
    std::shared_ptr<Transport_Discretization> transport_discretization;
    std::shared_ptr<Vector_Operator> source_operator;
    std::shared_ptr<Vector_Operator> flux_operator;
    std::shared_ptr<Vector_Operator> value_operator;
    
    
    // Output data
    Result result_;
};

#endif
