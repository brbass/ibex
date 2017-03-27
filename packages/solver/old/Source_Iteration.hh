#ifndef Source_Iteration_hh
#define Source_Iteration_hh

#include <memory>

#include "Solver.hh"
#include "Vector_Operator.hh"

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

/*
  Uses source iteration to solve the problem
  phi = D Linv M S phi + D Linv q

  D: Discrete to moment
  Linv: Inverse of transport operator
  M: Moment to discrete
  phi: Moment representation of the angular flux
  q: Discrete representation of the internal source
*/
class Source_Iteration : public Solver
{
public:

    // Constructor
    Source_Iteration(int max_iterations,
                     int solver_print,
                     double tolerance,
                     std::shared_ptr<Spatial_Discretization> spatial_discretization,
                     std::shared_ptr<Angular_Discretization> angular_discretization,
                     std::shared_ptr<Energy_Discretization> energy_discretization,
                     std::shared_ptr<Transport_Discretization> transport_discretization,
                     std::shared_ptr<Sweep_Operator> sweeper,
                     std::shared_ptr<Discrete_To_Moment> discrete_to_moment,
                     std::shared_ptr<Moment_To_Discrete> moment_to_discrete,
                     std::shared_ptr<Scattering> scattering,
                     std::shared_ptr<Fission> fission,
                     std::shared_ptr<Preconditioner> preconditioner = std::shared_ptr<Preconditioner>());
    
    // Solve fixed source problem
    virtual void solve_steady_state(vector<double> &x) override;

    // Solve k-eigenvalue problem
    virtual void solve_k_eigenvalue(double &k_eigenvalue, vector<double> &x) override;

    // Solve time-dependent problem
    virtual void solve_time_dependent(vector<double> &x) override;

    // Output data to XML file
    virtual void output(pugi::xml_node &output_node) const override;

private:

    // Check for convergence based on pointwise error in scalar flux
    bool check_phi_convergence(vector<double> const &x, 
                               vector<double> const &x_old,
                               double &error) const;
    
    bool check_k_convergence(double k,
                             double k_old,
                             double &error) const;
    
    double update_k(double k_eigenvalue_old,
                    vector<double> const &x,
                    vector<double> const &x_old) const;
    
    bool preconditioned_;
    int max_iterations_;
    int total_iterations_;
    int source_iterations_;
    double tolerance_;

    std::shared_ptr<Spatial_Discretization> spatial_discretization_;
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;
    std::shared_ptr<Transport_Discretization> transport_discretization_;
    std::shared_ptr<Sweep_Operator> sweeper_;
    std::shared_ptr<Discrete_To_Moment> discrete_to_moment_;
    std::shared_ptr<Moment_To_Discrete> moment_to_discrete_;
    std::shared_ptr<Scattering> scattering_;
    std::shared_ptr<Fission> fission_;
    std::shared_ptr<Preconditioner> preconditioner_;
    
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
        
        Source_Iterator(Source_Iteration const &source_iteration);

        virtual void check_class_invariants() const override
        {
        }
        
    private:
        
        virtual void apply(vector<double> &x) const override;
        
        Source_Iteration const &si_;
    };
    
    /*
      Computes the moments of the flux via source iteration
      phi = D Linv M S phi + b
      
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
        Flux_Iterator(Source_Iteration const &source_iteration);
        
        virtual void check_class_invariants() const override
        {
        }
        
    private:
        
        virtual void apply(vector<double> &x) const override;
        
        Source_Iteration const &si_;
    };
};


#endif
