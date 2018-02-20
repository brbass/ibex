#ifndef Meshless_Sweep_Parser_hh
#define Meshless_Sweep_Parser_hh

#include <memory>

class Angular_Discretization;
class Energy_Discretization;
class Transport_Discretization;
class Meshless_Sweep;
class Weak_Spatial_Discretization;
class XML_Node;

class Meshless_Sweep_Parser
{
public:
    
    Meshless_Sweep_Parser(std::shared_ptr<Weak_Spatial_Discretization> spatial,
                      std::shared_ptr<Angular_Discretization> angular,
                      std::shared_ptr<Energy_Discretization> energy,
                      std::shared_ptr<Transport_Discretization> transport);
    
    std::shared_ptr<Meshless_Sweep> get_meshless_sweep(XML_Node input_node) const;

private:
    
    std::shared_ptr<Meshless_Sweep> get_weak_sweep(XML_Node input_node) const;
    std::shared_ptr<Meshless_Sweep> get_strong_sweep(XML_Node input_node) const;
    
    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
    std::shared_ptr<Transport_Discretization> transport_;
};

#endif
