#ifndef Krylov_Iteration_hh
#define Krylov_Iteration_hh

#include <memory>

#include "Solver.hh"
#include "Vector_Operator.hh"

class Epetra_Comm;
class Epetra_Map;

class Angular_Discretization;
class Discrete_To_Moment;
class Energy_Discretization;
class Fission;
class Moment_To_Discrete;
class Preconditioner;
class Scattering;
class Spatial_Discretization;
class Sweep_Operator;
class Transport_Discretization;

using std::shared_ptr;

/*
  Uses Krylov methods to solve the problem 
  (I - D Linv M S) phi = D Linv q
  
  I: Identity
  D: Discrete to moment
  Linv: Inverse of transport operator
  M: Moment to discrete
  phi: Moment representation of the angular flux
  q: Discrete representation of the internal source
*/
class Krylov_Iteration : public Solver
{
public:
    
    // Constructor
    Krylov_Iteration(int max_iterations,
                     int kspace, // Number of past guesses to store
                     int solver_print, 
                     double tolerance,
                     shared_ptr<Spatial_Discretization> spatial_discretization,
                     shared_ptr<Angular_Discretization> angular_discretization,
                     shared_ptr<Energy_Discretization> energy_discretization,
                     shared_ptr<Transport_Discretization> transport_discretization,
                     shared_ptr<Sweep_Operator> sweeper,
                     shared_ptr<Discrete_To_Moment> discrete_to_moment,
                     shared_ptr<Moment_To_Discrete> moment_to_discrete,
                     shared_ptr<Scattering> scattering,
                     shared_ptr<Fission> fission,
                     shared_ptr<Preconditioner> preconditioner = shared_ptr<Preconditioner>());
    
    // Solve fixed source problem
    virtual void solve_steady_state(vector<double> &x) override;
    
    // Solve k-eigenvalue problem
    virtual void solve_k_eigenvalue(double &k_eigenvalue, vector<double> &x) override;
    
    // Solve time-dependent problem
    virtual void solve_time_dependent(vector<double> &x) override;
    
    // Output data to XML file
    virtual void output(pugi::xml_node &output_node) const override;
    
    // Check for convergence based on pointwise error in scalar flux
    bool check_phi_convergence(vector<double> const &x, 
                               vector<double> const &x_old,
                               double &error);
    
private:

    bool preconditioned_;
    int max_iterations_;
    int kspace_;
    int total_iterations_;
    int source_iterations_;
    double tolerance_;

    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
    shared_ptr<Transport_Discretization> transport_discretization_;
    shared_ptr<Sweep_Operator> sweeper_;
    shared_ptr<Discrete_To_Moment> discrete_to_moment_;
    shared_ptr<Moment_To_Discrete> moment_to_discrete_;
    shared_ptr<Scattering> scattering_;
    shared_ptr<Fission> fission_;
    shared_ptr<Preconditioner> preconditioner_;
    
    /*
      Computes the first-flight flux, including boundary sources, via source iteration
      b = D Linv q

      D: Discrete to moment
      Linv: Inverse of transport operator
      M: Moment to discrete
      q: Discrete representation of the internal source
      b: Moment representation of the first-flight flux
    */
    class Source_Iterator : public Vector_Operator
    {
    public:

        // Constructor
        Source_Iterator(Krylov_Iteration const &krylov_iteration);
        
        virtual void check_class_invariants() const override
        {
        }

    private:
        
        virtual void apply(vector<double> &x) const override;
        
        Krylov_Iteration const &ki_;
    };

    /*
      Iteratively computes the moments of the flux using Krylov methods
      (I - D Linv M S) phi = b

      I: Identity
      D: Discrete to moment
      Linv: Inverse of transport operator
      M: Moment to discrete
      phi: Moment representation of the angular flux
      b: Moment representation of the first-flight flux
    */
    class Flux_Iterator : public Vector_Operator
    {
    public:

        // Constructor
        Flux_Iterator(Krylov_Iteration const &krylov_iteration,
                      bool include_fission = true);
        
        virtual void check_class_invariants() const override
        {
        }

    private:
        
        virtual void apply(vector<double> &x) const override;

        bool include_fission_;
        Krylov_Iteration const &ki_;
    };

    /* 
       Represents the fission operator for an eigenvalue calculation,
       D Linv M F,
    */
    class Fission_Iterator : public Vector_Operator
    {
    public:

        // Constructor
        Fission_Iterator(Krylov_Iteration const &krylov_iteration);
        
        virtual void check_class_invariants() const override
        {
        }

    private:
        
        virtual void apply(vector<double> &x) const override;
        
        Krylov_Iteration const &ki_;
    };

    /* 
       Represents the eigenvalue operator,
       (I - D Linv M S)inv D Linv M F
    */
    class Eigenvalue_Iterator : public Vector_Operator
    {
    public:

        // Constructor
        Eigenvalue_Iterator(Krylov_Iteration const &krylov_iteration,
                            shared_ptr<Epetra_Comm> comm,
                            shared_ptr<Epetra_Map> map);
        
        virtual void check_class_invariants() const override
        {
        }

    private:
        
        virtual void apply(vector<double> &x) const override;
        
        Krylov_Iteration const &ki_;
        shared_ptr<Epetra_Comm> comm_;
        shared_ptr<Epetra_Map> map_;
    };
};


#endif
