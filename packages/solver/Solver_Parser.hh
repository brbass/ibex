#ifndef Solver_Parser_hh
#define Solver_Parser_hh

#include "Parser.hh"
#include "Solver.hh"

class Angular_Discretization;
class Discrete_To_Moment;
class Energy_Discretization;
class Fission;
class Krylov_Iteration;
class Moment_To_Discrete;
class Power_Iteration;
class Preconditioner;
class Sweep_Operator;
class Scattering;
class Scattering_Operator;
class Source_Iteration;
class Spatial_Discretization;
class Transport_Discretization;

/*
  Create a Solver object
*/
class Solver_Parser : public Parser<Solver>
{
public:

    // Constructor
    Solver_Parser(pugi::xml_node &input_file,
                  shared_ptr<Spatial_Discretization> spatial,
                  shared_ptr<Angular_Discretization> angular,
                  shared_ptr<Energy_Discretization> energy,
                  shared_ptr<Transport_Discretization> transport,
                  shared_ptr<Sweep_Operator> sweeper);

    // Return object
    virtual shared_ptr<Solver> get_ptr() override
    {
        return solver_;
    }

    // Parse source iteration
    shared_ptr<Source_Iteration> parse_source_iteration();

    // Parse Krylov iteration
    shared_ptr<Krylov_Iteration> parse_krylov_iteration();

    // Parse power iteration
    shared_ptr<Power_Iteration> parse_power_iteration();
    
    // Parse the Discrete_To_Moment operator
    shared_ptr<Discrete_To_Moment> parse_discrete_to_moment();

    // Parse the Moment_To_Discrete operator
    shared_ptr<Moment_To_Discrete> parse_moment_to_discrete();

    // Parse the scattering operator
    shared_ptr<Scattering> parse_scattering();

    // Parse the fission operator
    shared_ptr<Fission> parse_fission();
    
    // Parse the preconditioner
    shared_ptr<Preconditioner> parse_preconditioner(pugi::xml_node solver_node,
                                                    shared_ptr<Scattering_Operator> scattering,
                                                    shared_ptr<Scattering_Operator> fission);
    
private:
    
    shared_ptr<Solver> solver_;
    
    shared_ptr<Spatial_Discretization> spatial_;
    shared_ptr<Angular_Discretization> angular_;
    shared_ptr<Energy_Discretization> energy_;
    shared_ptr<Transport_Discretization> transport_;
    shared_ptr<Sweep_Operator> sweeper_;
};

#endif
