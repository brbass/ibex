#ifndef Krylov_Steady_State_hh
#define Krylov_Steady_State_hh

#include "Solver.hh"

#include <memory>

class Angular_Discretization;
class Convergence_Measure;
class Energy_Discretization;
class Spatial_Discretization;
class Transport_Discretization;
class Vector_Operator;

class Krylov_Steady_State : public Solver
{
public:

    struct Options
    {
        int max_source_iterations = 5000;
        int max_iterations = 1000;
        int kspace = 20; // Number of past guesses to store
        int solver_print = 0;
        double tolerance = 1e-10;
    };

    Krylov_Steady_State(Options options,
                        std::shared_ptr<Spatial_Discretization> spatial_discretization,
                        std::shared_ptr<Angular_Discretization> angular_discretization,
                        std::shared_ptr<Energy_Discretization> energy_discretization,
                        std::shared_ptr<Transport_Discretization> transport_discretization,
                        std::shared_ptr<Convergence_Measure> convergence,
                        std::shared_ptr<Vector_Operator> source_operator,
                        std::shared_ptr<Vector_Operator> flux_operator,
                        std::vector<std::shared_ptr<Vector_Operator> > value_operators);
    
    virtual void solve() override;
    virtual void output(XML_Node output_node) const override;
    virtual void check_class_invariants() const override;
    virtual std::shared_ptr<Result> result() const override
    {
        return result_;
    }
    
protected:
    
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
