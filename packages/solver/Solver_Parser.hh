#ifndef Solver_Parser_hh
#define Solver_Parser_hh

#include <memory>

class Angular_Discretization;
class Energy_Discretization;
class Krylov_Steady_State;
class Krylov_Eigenvalue;
class Solver;
class Solver_Factory;
class Source_Iteration;
class Sweep_Operator;
class Transport_Discretization;
class Weak_Spatial_Discretization;
class XML_Node;

class Solver_Parser
{
public:
    Solver_Parser(std::shared_ptr<Weak_Spatial_Discretization> spatial,
                  std::shared_ptr<Angular_Discretization> angular,
                  std::shared_ptr<Energy_Discretization> energy,
                  std::shared_ptr<Transport_Discretization> transport);

    std::shared_ptr<Source_Iteration>
    get_source_iteration(XML_Node input_node,
                         std::shared_ptr<Sweep_Operator> Linv) const;
    std::shared_ptr<Krylov_Steady_State>
    get_krylov_steady_state(XML_Node input_node,
                            std::shared_ptr<Sweep_Operator> Linv) const;
    std::shared_ptr<Krylov_Eigenvalue>
    get_krylov_eigenvalue(XML_Node input_node,
                          std::shared_ptr<Sweep_Operator> Linv) const;

private:
    
    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
    std::shared_ptr<Transport_Discretization> transport_;
    std::shared_ptr<Solver_Factory> factory_;
};

#endif
