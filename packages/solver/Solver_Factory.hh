#ifndef Solver_Factory_hh
#define Solver_Factory_hh

#include <memory>

class Angular_Discretization;
class Convergence_Measure;
class Energy_Discretization;
class Source_Iteration;
class Sweep_Operator;
class Transport_Discretization;
class Weak_Spatial_Discretization;

class Solver_Factory
{
public:
    
    Solver_Factory(std::shared_ptr<Weak_Spatial_Discretization> spatial,
                   std::shared_ptr<Angular_Discretization> angular,
                   std::shared_ptr<Energy_Discretization> energy,
                   std::shared_ptr<Transport_Discretization> transport);
    
    std::shared_ptr<Source_Iteration> get_supg_source_iteration(std::shared_ptr<Sweep_Operator> Linv,
                                                           std::shared_ptr<Convergence_Measure> convergence) const;
    
private:
    
    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
    std::shared_ptr<Transport_Discretization> transport_;
};

#endif
