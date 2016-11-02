#ifndef Sweep_Parser_hh
#define Sweep_Parser_hh

#include "Parser.hh"
#include "Sweep_Operator.hh"

class Angular_Discretization;
class Discrete_To_Moment;
class Energy_Discretization;
class Fission;
class Moment_To_Discrete;
class Preconditioner;
class Sweep_Operator;
class RBF_Collocation_Sweep;
class Scattering;
class Scattering_Operator;
class Spatial_Discretization;
class Transport_Discretization;

/*
  Create a Sweep object
*/
class Sweep_Parser : public Parser<Sweep_Operator>
{
public:
    
    // Constructor
    Sweep_Parser(pugi::xml_node &input_file,
                 shared_ptr<Spatial_Discretization> spatial,
                 shared_ptr<Angular_Discretization> angular,
                 shared_ptr<Energy_Discretization> energy,
                 shared_ptr<Transport_Discretization> transport);
    
    // Return object
    virtual shared_ptr<Sweep_Operator> get_ptr() override
    {
        return sweeper_;
    }
    
    // Parse the Sweep_Operator
    shared_ptr<Sweep_Operator> parse_sweeper();
    
    // Parse RBF_Collocation
    shared_ptr<RBF_Collocation_Sweep> parse_rbf_collocation(pugi::xml_node &input_node);
    
private:
    
    shared_ptr<Sweep_Operator> sweeper_;
    
    shared_ptr<Spatial_Discretization> spatial_;
    shared_ptr<Angular_Discretization> angular_;
    shared_ptr<Energy_Discretization> energy_;
    shared_ptr<Transport_Discretization> transport_;
};

#endif
