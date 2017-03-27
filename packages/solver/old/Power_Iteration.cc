#include "Power_Iteration.hh"

#include <cmath>
#include <iomanip>
#include <iostream>

#include <AztecOO.h>
#include <mpi.h>
#include <Epetra_MultiVector.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>

#include "Angular_Discretization.hh"
#include "Augmented_Operator.hh"
#include "Check.hh"
#include "Discrete_To_Moment.hh"
#include "Energy_Discretization.hh"
#include "Epetra_Operator_Interface.hh"
#include "Fission.hh"
#include "Identity_Operator.hh"
#include "Moment_To_Discrete.hh"
#include "Multiplicative_Operator.hh"
#include "Preconditioner.hh"
#include "Scattering.hh"
#include "Spatial_Discretization.hh"
#include "Sweep_Operator.hh"
#include "Transport_Discretization.hh"
#include "Vector_Functions.hh"
#include "Vector_Operator_Functions.hh"
#include "XML_Functions.hh"

using namespace std;

namespace vf = Vector_Functions;

Power_Iteration::
Power_Iteration(int max_iterations,
                int kspace,
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
                shared_ptr<Preconditioner> preconditioner):
    Solver(solver_print),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    transport_discretization_(transport_discretization),
    max_iterations_(max_iterations),
    kspace_(kspace),
    total_iterations_(0),
    tolerance_(tolerance),
    sweeper_(sweeper),
    discrete_to_moment_(discrete_to_moment),
    moment_to_discrete_(moment_to_discrete),
    scattering_(scattering),
    fission_(fission),
    preconditioner_(preconditioner)
{
    if (preconditioner_ == shared_ptr<Preconditioner>())
    {
        preconditioned_ = false;
    }
    else
    {
        preconditioned_ = true;
    }
}

void Power_Iteration::
solve_steady_state(vector<double> &x)
{
    AssertMsg(false, "not implemented");
}

void Power_Iteration::
solve_k_eigenvalue(double &k_eigenvalue, 
                   vector<double> &x)
{
    int phi_size = transport_discretization_->phi_size();
    int number_of_augments = transport_discretization_->number_of_augments();

    k_eigenvalue = 1.0;
    x.resize(phi_size + number_of_augments, 1.0);

    shared_ptr<Epetra_Comm> comm = make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    shared_ptr<Epetra_Map> map = make_shared<Epetra_Map>(phi_size + number_of_augments, 0, *comm);
    
    shared_ptr<Vector_Operator> FI
        = make_shared<Flux_Iterator>(*this,
                                     false); // don't include fission
    shared_ptr<Vector_Operator> NI
        = make_shared<Fission_Iterator>(*this);
    shared_ptr<Vector_Operator> EI
        = make_shared<Eigenvalue_Iterator>(*this,
                                           comm,
                                           map);
    shared_ptr<Multiplicative_Operator> K
        = make_shared<Multiplicative_Operator>(phi_size + number_of_augments,
                                               k_eigenvalue);
    
    shared_ptr<Vector_Operator> Op = EI * K;
    
    double k_eigenvalue_old;
    vector<double> x_old;
    
    print_name("Power iteration");
    
    for (int it = 0; it < max_iterations_; ++it)
    {
        print_iteration(it);
        
        k_eigenvalue_old = k_eigenvalue;
        x_old = x;
        
        K->set_scalar(1. / k_eigenvalue);
        
        (*Op)(x);
        
        k_eigenvalue = update_k(k_eigenvalue_old,
                                x,
                                x_old);
        double error;
        bool converged = check_k_convergence(k_eigenvalue, k_eigenvalue_old, error);
        print_value(k_eigenvalue);
        print_error(error);
        
        if (converged)
        {
            total_iterations_ = it + 1;
            
            print_convergence();
            print_eigenvalue(k_eigenvalue);
            
            break;
        }
    }
    
    if (total_iterations_ == max_iterations_)
    {
        print_failure();
    }
    
    x.resize(phi_size); // remove augments
}

void Power_Iteration::
solve_time_dependent(vector<double> &x)
{
    AssertMsg(false, "not implemented");
}

void Power_Iteration::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node source = output_node.append_child("solution_method");

    XML_Functions::append_child(source, "power_iteration", "type");
    XML_Functions::append_child(source, max_iterations_, "max_iterations");
    XML_Functions::append_child(source, total_iterations_, "total_iterations");
    XML_Functions::append_child(source, tolerance_, "tolerance");
}

double Power_Iteration::
update_k(double k_eigenvalue_old,
         vector<double> const &x,
         vector<double> const &x_old) const
{
    int phi_size = transport_discretization_->phi_size();

    vector<double> fis(x.begin(), x.begin() + phi_size);
    vector<double> fis_old(x_old.begin(), x_old.begin() + phi_size);
    
    (*fission_)(fis);
    (*fission_)(fis_old);
    
    double num = vf::magnitude(fis) * k_eigenvalue_old;
    double den = vf::magnitude(fis_old);
    
    return num / den;
}

bool Power_Iteration::
check_k_convergence(double k,
                    double k_old,
                    double &error) const
{
    error = abs(k - k_old) / (abs(k_old) + tolerance_ * tolerance_);
    
    if (error > tolerance_)
    {
        return false;
    }
    else
    {
        return true;
    }
}

Power_Iteration::Flux_Iterator::
Flux_Iterator(Power_Iteration const &pi,
              bool include_fission):
    Vector_Operator(pi.transport_discretization_->phi_size() + pi.transport_discretization_->number_of_augments(),
                    pi.transport_discretization_->phi_size() + pi.transport_discretization_->number_of_augments()),
    include_fission_(include_fission),
    pi_(pi)
{
    
}

void Power_Iteration::Flux_Iterator::
apply(vector<double> &x) const
{
    int phi_size = pi_.transport_discretization_->phi_size();
    int number_of_augments = pi_.transport_discretization_->number_of_augments();
    
    shared_ptr<Vector_Operator> D
        = make_shared<Augmented_Operator>(number_of_augments,
                                          pi_.discrete_to_moment_);
    shared_ptr<Vector_Operator> M
        = make_shared<Augmented_Operator>(number_of_augments,
                                          pi_.moment_to_discrete_);

    shared_ptr<Sweep_Operator> Linv = pi_.sweeper_;
    Assert(Linv);
    Linv->set_include_boundary_source(false);

    shared_ptr<Vector_Operator> T;
    switch (Linv->sweep_type())
    {
    case Sweep_Operator::Sweep_Type::MOMENT:
        T = Linv;
        break;
    case Sweep_Operator::Sweep_Type::ORDINATE:
        T = D * Linv * M;
        break;
    }
    
    shared_ptr<Vector_Operator> E;
    if (pi_.preconditioned_)
    {
        shared_ptr<Preconditioner> Cinv  = pi_.preconditioner_;
        Assert(Cinv);
        shared_ptr<Sweep_Operator> sweeper = Cinv->sweeper();
        sweeper->set_include_boundary_source(false);
        switch (sweeper->sweep_type())
        {
        case Sweep_Operator::Sweep_Type::MOMENT:
            E = Cinv;
            break;
        case Sweep_Operator::Sweep_Type::ORDINATE:
            E = D * Cinv * M;
            break;
        }
    }
    
    // Add augments to operators
    
    shared_ptr<Vector_Operator> S
        = make_shared<Augmented_Operator>(number_of_augments,
                                          pi_.scattering_);
    shared_ptr<Vector_Operator> SF;
    
    if (include_fission_)
    {
        shared_ptr<Vector_Operator> F
            = make_shared<Augmented_Operator>(number_of_augments,
                                              pi_.fission_,
                                              true); // Only one of S or F needs to include the augments
        
        SF = S + F;
    }
    else
    {
        SF = S;
    }

    shared_ptr<Vector_Operator> I
        = make_shared<Identity_Operator>(phi_size + number_of_augments);
    
    // Create combined operators
    
    shared_ptr<Vector_Operator> Op;

    if (pi_.preconditioned_)
    {
        Op = (I + E * SF) * (I - T * SF);
    }
    else
    {
        Op = I - T * SF;
    }
    
    (*Op)(x);
}

Power_Iteration::Fission_Iterator::
Fission_Iterator(Power_Iteration const &pi):
    Vector_Operator(pi.transport_discretization_->phi_size() + pi.transport_discretization_->number_of_augments(),
                    pi.transport_discretization_->phi_size() + pi.transport_discretization_->number_of_augments()),
    pi_(pi)
{
}

void Power_Iteration::Fission_Iterator::
apply(vector<double> &x) const
{
    int number_of_augments = pi_.transport_discretization_->number_of_augments();
    
    shared_ptr<Vector_Operator> D
        = make_shared<Augmented_Operator>(number_of_augments,
                                          pi_.discrete_to_moment_);
    shared_ptr<Vector_Operator> M
        = make_shared<Augmented_Operator>(number_of_augments,
                                          pi_.moment_to_discrete_);
    
    shared_ptr<Sweep_Operator> Linv = pi_.sweeper_;
    Linv->set_include_boundary_source(false);
    
    shared_ptr<Vector_Operator> T;
    switch (Linv->sweep_type())
    {
    case Sweep_Operator::Sweep_Type::MOMENT:
        T = Linv;
        break;
    case Sweep_Operator::Sweep_Type::ORDINATE:
        T = D * Linv * M;
        break;
    }
    
    shared_ptr<Vector_Operator> F
        = make_shared<Augmented_Operator>(number_of_augments,
                                          pi_.fission_,
                                          true); // zero out augments for fission source
    
    
    shared_ptr<Vector_Operator> Op
        = T * F;

    (*Op)(x);
}

Power_Iteration::Eigenvalue_Iterator::
Eigenvalue_Iterator(Power_Iteration const &pi,
                    shared_ptr<Epetra_Comm> comm,
                    shared_ptr<Epetra_Map> map):
    Vector_Operator(pi.transport_discretization_->phi_size() + pi.transport_discretization_->number_of_augments(),
                    pi.transport_discretization_->phi_size() + pi.transport_discretization_->number_of_augments()),
    pi_(pi),
    comm_(comm),
    map_(map)
{
}

void Power_Iteration::Eigenvalue_Iterator::
apply(vector<double> &x) const
{
    shared_ptr<Vector_Operator> FI
        = make_shared<Flux_Iterator>(pi_,
                                     false); // don't include fission
    shared_ptr<Vector_Operator> NI
        = make_shared<Fission_Iterator>(pi_);
    
    (*NI)(x);
    
    shared_ptr<Epetra_Vector> lhs
        = make_shared<Epetra_Vector>(*map_);
    shared_ptr<Epetra_Vector> rhs
        = make_shared<Epetra_Vector>(Copy,
                                     *map_,
                                     &x[0]);
    shared_ptr<Epetra_Operator> oper
        = make_shared<Epetra_Operator_Interface>(comm_,
                                                 map_,
                                                 FI);
    shared_ptr<Epetra_LinearProblem> problem
        = make_shared<Epetra_LinearProblem>(oper.get(),
                                            lhs.get(),
                                            rhs.get());
    shared_ptr<AztecOO> solver
        = make_shared<AztecOO>(*problem);
    
    solver->SetAztecOption(AZ_precond, AZ_none);
    solver->SetAztecOption(AZ_solver, AZ_gmres);
    solver->SetAztecOption(AZ_kspace, pi_.kspace_);
    solver->SetAztecOption(AZ_conv, AZ_rhs);
    if (pi_.solver_print_ && false)
    {
        solver->SetAztecOption(AZ_output, AZ_all);
    }
    else
    {
        solver->SetAztecOption(AZ_output, AZ_none);
    }
    
    lhs->PutScalar(1.0);
    
    solver->Iterate(pi_.max_iterations_, pi_.tolerance_);
    
    lhs->ExtractCopy(&x[0]);
}

