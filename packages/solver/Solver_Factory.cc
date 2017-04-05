#include "Solver_Factory.hh"

#include "Angular_Discretization.hh"
#include "Augmented_Operator.hh"
#include "Boundary_Source_Toggle.hh"
#include "Discrete_To_Moment.hh"
#include "Discrete_Normalization_Operator.hh"
#include "Energy_Discretization.hh"
#include "Fission.hh"
#include "Identity_Operator.hh"
#include "Internal_Source_Operator.hh"
#include "Krylov_Eigenvalue.hh"
#include "Krylov_Steady_State.hh"
#include "Moment_To_Discrete.hh"
#include "Moment_Value_Operator.hh"
#include "Moment_Weighting_Operator.hh"
#include "Resize_Operator.hh"
#include "Scattering.hh"
#include "Source_Iteration.hh"
#include "SUPG_Fission.hh"
#include "SUPG_Internal_Source_Operator.hh"
#include "SUPG_Moment_To_Discrete.hh"
#include "SUPG_Scattering.hh"
#include "Transport_Discretization.hh"
#include "Vector_Operator_Functions.hh"
#include "Weak_Spatial_Discretization.hh"

using std::make_shared;
using std::shared_ptr;
using std::vector;

Solver_Factory::
Solver_Factory(shared_ptr<Weak_Spatial_Discretization> spatial,
               shared_ptr<Angular_Discretization> angular,
               shared_ptr<Energy_Discretization> energy,
               shared_ptr<Transport_Discretization> transport):
    spatial_(spatial),
    angular_(angular),
    energy_(energy),
    transport_(transport)
{
}

void Solver_Factory::
get_source_operators(shared_ptr<Sweep_Operator> Linv,
                     shared_ptr<Vector_Operator> &source_operator,
                     shared_ptr<Vector_Operator> &flux_operator) const
{
    // Check if problem includes SUPG terms
    bool include_supg = spatial_->include_supg();
    if (include_supg)
    {
        return get_supg_source_operators(Linv,
                                         source_operator,
                                         flux_operator);
                                         
    }

    // Get size data
    int phi_size = transport_->phi_size();
    int number_of_augments = transport_->number_of_augments();
    
    // Get moment-to-discrete and discrete-to-moment operators
    shared_ptr<Vector_Operator> M
        = make_shared<Moment_To_Discrete>(spatial_,
                                          angular_,
                                          energy_);
    shared_ptr<Vector_Operator> D
        = make_shared<Discrete_To_Moment>(spatial_,
                                          angular_,
                                          energy_);
    
    // Get Scattering operators
    Scattering_Operator::Options scattering_options;
    shared_ptr<Vector_Operator> S
        = make_shared<Scattering>(spatial_,
                                  angular_,
                                  energy_,
                                  scattering_options);
    shared_ptr<Vector_Operator> F
        = make_shared<Fission>(spatial_,
                               angular_,
                               energy_,
                               scattering_options);
    
    // Get weighting operator
    Weighting_Operator::Options weighting_options;
    shared_ptr<Vector_Operator> Wm
        = make_shared<Moment_Weighting_Operator>(spatial_,
                                                 angular_,
                                                 energy_,
                                                 weighting_options);
    
    // Get source operator
    shared_ptr<Vector_Operator> Q
        = make_shared<Internal_Source_Operator>(spatial_,
                                                angular_,
                                                energy_);
    
    // Add augments to operators
    if (number_of_augments > 0)
    {
        M = make_shared<Augmented_Operator>(number_of_augments,
                                            M,
                                            false);
        D = make_shared<Augmented_Operator>(number_of_augments,
                                            D,
                                            false);
        S = make_shared<Augmented_Operator>(number_of_augments,
                                            S,
                                            false);
        F = make_shared<Augmented_Operator>(number_of_augments,
                                            F,
                                            true);
        Wm = make_shared<Augmented_Operator>(number_of_augments,
                                             Wm,
                                             false);
        Q = make_shared<Augmented_Operator>(number_of_augments,
                                            Q,
                                            false);
    }
    
    // Get sweep operator with boundary source off/on
    shared_ptr<Vector_Operator> LinvB
        = make_shared<Boundary_Source_Toggle>(true,
                                              Linv);
    shared_ptr<Vector_Operator> LinvI
        = make_shared<Boundary_Source_Toggle>(false,
                                              Linv);
    
    // Get combined operators
    source_operator
        = D * LinvB * M * Q;
    flux_operator
        = D * LinvI * M * (S + F) * Wm;
}

void Solver_Factory::
get_supg_source_operators(shared_ptr<Sweep_Operator> Linv,
                          shared_ptr<Vector_Operator> &source_operator,
                          shared_ptr<Vector_Operator> &flux_operator) const
{
    // Check that this problem is SUPG
    bool include_supg = spatial_->include_supg();
    Assert(include_supg);

    // Get size data
    int number_of_dimensional_moments = spatial_->number_of_dimensional_moments();
    int phi_size = transport_->phi_size();
    int number_of_augments = transport_->number_of_augments();
    
    // Get moment-to-discrete and discrete-to-moment operators
    shared_ptr<Vector_Operator> M1
        = make_shared<SUPG_Moment_To_Discrete>(spatial_,
                                               angular_,
                                               energy_,
                                               false); // include double dimenisonal moments
    shared_ptr<Vector_Operator> M2
        = make_shared<SUPG_Moment_To_Discrete>(spatial_,
                                               angular_,
                                               energy_,
                                               true); // include double dimensional moments
    shared_ptr<Vector_Operator> D
        = make_shared<Discrete_To_Moment>(spatial_,
                                          angular_,
                                          energy_);
    
    // Get Scattering operators
    SUPG_Scattering_Operator::Options scattering_options;
    shared_ptr<Vector_Operator> S
        = make_shared<SUPG_Scattering>(spatial_,
                                       angular_,
                                       energy_,
                                       scattering_options);
    shared_ptr<Vector_Operator> F
        = make_shared<SUPG_Fission>(spatial_,
                                    angular_,
                                    energy_,
                                    scattering_options);
    
    // Get normalization operator
    Weighting_Operator::Options weighting_options;
    shared_ptr<Vector_Operator> W
        = make_shared<Moment_Weighting_Operator>(spatial_,
                                                 angular_,
                                                 energy_,
                                                 weighting_options);
    shared_ptr<Vector_Operator> N
        = make_shared<Discrete_Normalization_Operator>(spatial_,
                                                       angular_,
                                                       energy_,
                                                       weighting_options);
    
    // Get source operator
    shared_ptr<Vector_Operator> Q
        = make_shared<SUPG_Internal_Source_Operator>(spatial_,
                                                     angular_,
                                                     energy_);
    
    // Add augments to operators
    if (number_of_augments > 0)
    {
        M1 = make_shared<Augmented_Operator>(number_of_augments,
                                             M1,
                                             false);
        M2 = make_shared<Augmented_Operator>(number_of_augments,
                                             M2,
                                             false);
        D = make_shared<Augmented_Operator>(number_of_augments,
                                            D,
                                            false);
        S = make_shared<Augmented_Operator>(number_of_augments,
                                            S,
                                            false);
        F = make_shared<Augmented_Operator>(number_of_augments,
                                            F,
                                            true);
        W = make_shared<Augmented_Operator>(number_of_augments,
                                            W,
                                            false);
        N = make_shared<Augmented_Operator>(number_of_augments,
                                            N,
                                            false);
        Q = make_shared<Augmented_Operator>(number_of_augments,
                                            Q,
                                            false);
    }
    
    // Get sweep operator with boundary source off/on
    shared_ptr<Vector_Operator> LinvB
        = make_shared<Boundary_Source_Toggle>(true,
                                              Linv);
    shared_ptr<Vector_Operator> LinvI
        = make_shared<Boundary_Source_Toggle>(false,
                                              Linv);
    
    // Get combined operators
    source_operator
        = D * LinvB * M1 * Q;
    flux_operator
        = D * LinvI * N * M2 * (S + F) * W;
}

void Solver_Factory::
get_eigenvalue_operators(shared_ptr<Sweep_Operator> Linv,
                         shared_ptr<Vector_Operator> &fission_operator,
                         shared_ptr<Vector_Operator> &flux_operator) const
{
    // Check if problem includes SUPG terms
    bool include_supg = spatial_->include_supg();
    if (include_supg)
    {
        return get_supg_eigenvalue_operators(Linv,
                                             fission_operator,
                                             flux_operator);
                                         
    }

    // Get size data
    int phi_size = transport_->phi_size();
    int number_of_augments = transport_->number_of_augments();
    
    // Get moment-to-discrete and discrete-to-moment operators
    shared_ptr<Vector_Operator> M
        = make_shared<Moment_To_Discrete>(spatial_,
                                          angular_,
                                          energy_);
    shared_ptr<Vector_Operator> D
        = make_shared<Discrete_To_Moment>(spatial_,
                                          angular_,
                                          energy_);
    
    // Get Scattering operators
    Scattering_Operator::Options scattering_options;
    shared_ptr<Vector_Operator> S
        = make_shared<Scattering>(spatial_,
                                  angular_,
                                  energy_,
                                  scattering_options);
    shared_ptr<Vector_Operator> F
        = make_shared<Fission>(spatial_,
                               angular_,
                               energy_,
                               scattering_options);
    
    // Get weighting operator
    Weighting_Operator::Options weighting_options;
    shared_ptr<Vector_Operator> W
        = make_shared<Moment_Weighting_Operator>(spatial_,
                                                 angular_,
                                                 energy_,
                                                 weighting_options);
    
    // Add augments to operators
    if (number_of_augments > 0)
    {
        M = make_shared<Augmented_Operator>(number_of_augments,
                                            M,
                                            false);
        D = make_shared<Augmented_Operator>(number_of_augments,
                                            D,
                                            false);
        S = make_shared<Augmented_Operator>(number_of_augments,
                                            S,
                                            false);
        F = make_shared<Augmented_Operator>(number_of_augments,
                                            F,
                                            true);
        W = make_shared<Augmented_Operator>(number_of_augments,
                                            W,
                                            false);
    }
    
    // Get sweep operator with boundary source off
    shared_ptr<Vector_Operator> LinvI
        = make_shared<Boundary_Source_Toggle>(false,
                                              Linv);
    
    // Get combined operators
    fission_operator
        = D * LinvI * M * F * W;
    flux_operator
        = D * LinvI * M * S * W;
}

void Solver_Factory::
get_supg_eigenvalue_operators(shared_ptr<Sweep_Operator> Linv,
                              shared_ptr<Vector_Operator> &fission_operator,
                              shared_ptr<Vector_Operator> &flux_operator) const
{
    // Check that this problem is SUPG
    bool include_supg = spatial_->include_supg();
    Assert(include_supg);

    // Get size data
    int number_of_dimensional_moments = spatial_->number_of_dimensional_moments();
    int phi_size = transport_->phi_size();
    int number_of_augments = transport_->number_of_augments();
    
    // Get moment-to-discrete and discrete-to-moment operators
    shared_ptr<Vector_Operator> M2
        = make_shared<SUPG_Moment_To_Discrete>(spatial_,
                                               angular_,
                                               energy_,
                                               true); // include double dimensional moments
    shared_ptr<Vector_Operator> D
        = make_shared<Discrete_To_Moment>(spatial_,
                                          angular_,
                                          energy_);
    
    // Get Scattering operators
    SUPG_Scattering_Operator::Options scattering_options;
    shared_ptr<Vector_Operator> S
        = make_shared<SUPG_Scattering>(spatial_,
                                       angular_,
                                       energy_,
                                       scattering_options);
    shared_ptr<Vector_Operator> F
        = make_shared<SUPG_Fission>(spatial_,
                                    angular_,
                                    energy_,
                                    scattering_options);
    
    // Get normalization operator
    Weighting_Operator::Options weighting_options;
    shared_ptr<Vector_Operator> W
        = make_shared<Moment_Weighting_Operator>(spatial_,
                                                 angular_,
                                                 energy_,
                                                 weighting_options);
    shared_ptr<Vector_Operator> N
        = make_shared<Discrete_Normalization_Operator>(spatial_,
                                                       angular_,
                                                       energy_,
                                                       weighting_options);
    
    // Add augments to operators
    if (number_of_augments > 0)
    {
        M2 = make_shared<Augmented_Operator>(number_of_augments,
                                             M2,
                                             false);
        D = make_shared<Augmented_Operator>(number_of_augments,
                                            D,
                                            false);
        S = make_shared<Augmented_Operator>(number_of_augments,
                                            S,
                                            false);
        F = make_shared<Augmented_Operator>(number_of_augments,
                                            F,
                                            true);
        W = make_shared<Augmented_Operator>(number_of_augments,
                                            W,
                                            false);
        N = make_shared<Augmented_Operator>(number_of_augments,
                                            N,
                                            false);
    }
    
    // Get sweep operator with boundary source off
    shared_ptr<Vector_Operator> LinvI
        = make_shared<Boundary_Source_Toggle>(false,
                                              Linv);
    
    // Get combined operators
    fission_operator
        = D * LinvI * N * M2 * F * W;
    flux_operator
        = D * LinvI * N * M2 * S * W;
}

std::shared_ptr<Source_Iteration> Solver_Factory::
get_source_iteration(shared_ptr<Sweep_Operator> Linv,
                     shared_ptr<Convergence_Measure> convergence) const
{
    // Get combined operators
    shared_ptr<Vector_Operator> source_operator;
    shared_ptr<Vector_Operator> flux_operator;
    get_source_operators(Linv,
                         source_operator,
                         flux_operator);
    
    // Get value operators
    vector<shared_ptr<Vector_Operator> > value_operators
        = {make_shared<Moment_Value_Operator>(spatial_,
                                              angular_,
                                              energy_,
                                              false), // no weighting
           make_shared<Moment_Value_Operator>(spatial_,
                                              angular_,
                                              energy_,
                                              true)}; // weighting
    
    // Get source iteration
    Source_Iteration::Options iteration_options;
    iteration_options.solver_print = 0;
    return make_shared<Source_Iteration>(iteration_options,
                                         spatial_,
                                         angular_,
                                         energy_,
                                         transport_,
                                         convergence,
                                         source_operator,
                                         flux_operator,
                                         value_operators); 
    
}

std::shared_ptr<Krylov_Steady_State> Solver_Factory::
get_krylov_steady_state(shared_ptr<Sweep_Operator> Linv,
                        shared_ptr<Convergence_Measure> convergence) const
{
    // Get combined operators
    shared_ptr<Vector_Operator> source_operator;
    shared_ptr<Vector_Operator> flux_operator;
    get_source_operators(Linv,
                         source_operator,
                         flux_operator);
    shared_ptr<Identity_Operator> identity
        = make_shared<Identity_Operator>(flux_operator->column_size());
    flux_operator = identity - flux_operator;
    
    // Get value operators
    vector<shared_ptr<Vector_Operator> > value_operators
        = {make_shared<Moment_Value_Operator>(spatial_,
                                              angular_,
                                              energy_,
                                              false), // no weighting
           make_shared<Moment_Value_Operator>(spatial_,
                                              angular_,
                                              energy_,
                                              true)}; // weighting
    
    // Get source iteration
    Krylov_Steady_State::Options iteration_options;
    iteration_options.solver_print = 0;
    return make_shared<Krylov_Steady_State>(iteration_options,
                                            spatial_,
                                            angular_,
                                            energy_,
                                            transport_,
                                            convergence,
                                            source_operator,
                                            flux_operator,
                                            value_operators); 
    
}

std::shared_ptr<Krylov_Eigenvalue> Solver_Factory::
get_krylov_eigenvalue(shared_ptr<Sweep_Operator> Linv) const
{
    // Get combined operators
    shared_ptr<Vector_Operator> fission_operator;
    shared_ptr<Vector_Operator> flux_operator;
    get_eigenvalue_operators(Linv,
                             fission_operator,
                             flux_operator);
    shared_ptr<Identity_Operator> identity
        = make_shared<Identity_Operator>(flux_operator->column_size());
    flux_operator = identity - flux_operator;
    
    // Get value operators
    vector<shared_ptr<Vector_Operator> > value_operators
        = {make_shared<Moment_Value_Operator>(spatial_,
                                              angular_,
                                              energy_,
                                              false), // no weighting
           make_shared<Moment_Value_Operator>(spatial_,
                                              angular_,
                                              energy_,
                                              true)}; // weighting
    
    // Get source iteration
    Krylov_Eigenvalue::Options iteration_options;
    iteration_options.solver_print = 0;
    return make_shared<Krylov_Eigenvalue>(iteration_options,
                                          spatial_,
                                          angular_,
                                          energy_,
                                          transport_,
                                          fission_operator,
                                          flux_operator,
                                          value_operators); 
}
