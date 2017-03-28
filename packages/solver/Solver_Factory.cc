#include "Solver_Factory.hh"

#include "Angular_Discretization.hh"
#include "Augmented_Operator.hh"
#include "Boundary_Source_Toggle.hh"
#include "Dimensional_Moment_Copy.hh"
#include "Dimensional_Moment_Summation.hh"
#include "Discrete_To_Moment.hh"
#include "Discrete_Weighting_Operator.hh"
#include "Energy_Discretization.hh"
#include "Fission.hh"
#include "Internal_Source_Operator.hh"
#include "Moment_To_Discrete.hh"
#include "Moment_Value_Operator.hh"
#include "Moment_Weighting_Operator.hh"
#include "Resize_Operator.hh"
#include "Scattering.hh"
#include "Source_Iteration.hh"
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

std::shared_ptr<Source_Iteration> Solver_Factory::
get_supg_source_iteration(shared_ptr<Sweep_Operator> Linv,
                          shared_ptr<Convergence_Measure> convergence) const
{
    // Check to be sure that this problem is SUPG
    int number_of_dimensional_moments = spatial_->number_of_dimensional_moments();
    int phi_size = transport_->phi_size();
    int number_of_augments = transport_->number_of_augments();
    Assert(number_of_dimensional_moments > 1);
    
    // Get moment-to-discrete and discrete-to-moment operators
    shared_ptr<Vector_Operator> M
        = make_shared<Moment_To_Discrete>(spatial_,
                                          angular_,
                                          energy_,
                                          true); // include dimensional moments
    shared_ptr<Vector_Operator> D
        = make_shared<Discrete_To_Moment>(spatial_,
                                          angular_,
                                          energy_,
                                          false); // include dimensional moments

    // Get Scattering operators
    Scattering_Operator::Options scattering_options;
    scattering_options.include_dimensional_moments = true;
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
    shared_ptr<Vector_Operator> Wd
        = make_shared<Discrete_Weighting_Operator>(spatial_,
                                                   angular_,
                                                   energy_,
                                                   weighting_options);
    
    // Get dimensional summation and copy operators
    shared_ptr<Vector_Operator> Ns
        = make_shared<Dimensional_Moment_Summation>(spatial_,
                                                    angular_,
                                                    energy_);
    shared_ptr<Vector_Operator> Nc
        = make_shared<Dimensional_Moment_Copy>(spatial_,
                                               angular_,
                                               energy_);

    // Get source operator
    shared_ptr<Vector_Operator> Q
        = make_shared<Internal_Source_Operator>(spatial_,
                                                angular_,
                                                energy_);

    // Get the operator to resize given source to correct size
    shared_ptr<Vector_Operator> R
        = make_shared<Resize_Operator>(phi_size * number_of_dimensional_moments + number_of_augments,
                                       phi_size + number_of_augments);

    
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
        Wd = make_shared<Augmented_Operator>(number_of_augments,
                                             Wd,
                                             false);
        Ns = make_shared<Augmented_Operator>(number_of_augments,
                                             Ns,
                                             false);
        Nc = make_shared<Augmented_Operator>(number_of_augments,
                                             Nc,
                                             false);
        Q = make_shared<Augmented_Operator>(number_of_augments,
                                            Q,
                                            false);
        R = make_shared<Augmented_Operator>(number_of_augments,
                                            R,
                                            false);
    }
    
    // Get sweep operator with boundary source off/on
    shared_ptr<Vector_Operator> LinvB
        = make_shared<Boundary_Source_Toggle>(true,
                                              Linv);
    shared_ptr<Vector_Operator> LinvI
        = make_shared<Boundary_Source_Toggle>(false,
                                              Linv);

    // Get value operators
    vector<shared_ptr<Vector_Operator> > value_operators
        = {make_shared<Moment_Value_Operator>(spatial_,
                                              angular_,
                                              energy_,
                                              false),
           make_shared<Moment_Value_Operator>(spatial_,
                                              angular_,
                                              energy_,
                                              true)};

    // Get combined operators
    shared_ptr<Vector_Operator> source_operator
        = D * LinvB * Ns * M * Q * R;
    shared_ptr<Vector_Operator> flux_operator
        = D * LinvI * Wd * Ns * M * (S + F) * Nc;
    
    // Get source iteration
    Source_Iteration::Options iteration_options;
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

