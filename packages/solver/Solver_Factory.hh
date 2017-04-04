#ifndef Solver_Factory_hh
#define Solver_Factory_hh

#include <memory>

class Angular_Discretization;
class Convergence_Measure;
class Energy_Discretization;
class Krylov_Steady_State;
class Source_Iteration;
class Sweep_Operator;
class Transport_Discretization;
class Vector_Operator;
class Weak_Spatial_Discretization;

class Solver_Factory
{
public:
    
    Solver_Factory(std::shared_ptr<Weak_Spatial_Discretization> spatial,
                   std::shared_ptr<Angular_Discretization> angular,
                   std::shared_ptr<Energy_Discretization> energy,
                   std::shared_ptr<Transport_Discretization> transport);

    void get_source_operators(std::shared_ptr<Sweep_Operator> Linv,
                              std::shared_ptr<Vector_Operator> &source_operator,
                              std::shared_ptr<Vector_Operator> &flux_operator) const;
    
    void get_supg_source_operators(std::shared_ptr<Sweep_Operator> Linv,
                                       std::shared_ptr<Vector_Operator> &source_operator,
                                       std::shared_ptr<Vector_Operator> &flux_operator) const;
    
    std::shared_ptr<Source_Iteration> get_source_iteration(std::shared_ptr<Sweep_Operator> Linv,
                                                           std::shared_ptr<Convergence_Measure> convergence) const;
    std::shared_ptr<Krylov_Steady_State> get_krylov_steady_state(std::shared_ptr<Sweep_Operator> Linv,
                                                                 std::shared_ptr<Convergence_Measure> convergence) const;
private:
    
    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
    std::shared_ptr<Transport_Discretization> transport_;
};

#endif
