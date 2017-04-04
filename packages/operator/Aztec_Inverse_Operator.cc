#include "Aztec_Inverse_Operator.hh"

#include <AztecOO.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>

#include "Epetra_Operator_Interface.hh"

using std::make_shared;
using std::shared_ptr;
using std::string;
using std::vector;

Aztec_Inverse_Operator::
Aztec_Inverse_Operator(Options options,
                       shared_ptr<Vector_Operator> vector_operator):
    Inverse_Operator(vector_operator),
    options_(options)
{
    initialize_trilinos();
    check_class_invariants();
}

void Aztec_Inverse_Operator::
initialize_trilinos()
{
    comm_ = make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    map_ = make_shared<Epetra_Map>(size_, 0, *comm_);
    lhs_ = make_shared<Epetra_Vector>(*map_);
    rhs_ = make_shared<Epetra_Vector>(*map_);
    oper_ = make_shared<Epetra_Operator_Interface>(comm_,
                                                   map_,
                                                   vector_operator_);
    problem_ = make_shared<Epetra_LinearProblem>(oper_.get(),
                                                 lhs_.get(),
                                                 rhs_.get());
    solver_ = make_shared<AztecOO>(*problem_);
    solver_->SetAztecOption(AZ_precond, AZ_none);
    solver_->SetAztecOption(AZ_solver, AZ_gmres);
    solver_->SetAztecOption(AZ_kspace, options_.kspace);
    solver_->SetAztecOption(AZ_conv, AZ_rhs);
    if (options_.solver_print)
    {
        solver_->SetAztecOption(AZ_output, AZ_all);
    }
    else
    {
        solver_->SetAztecOption(AZ_output, AZ_none);
    }
}

void Aztec_Inverse_Operator::
apply(vector<double> &x) const
{
    // Set rhs vector
    // Use x as initial guess
    for (int i = 0; i < size_; ++i)
    {
        (*lhs_)[i] = x[i];
        (*rhs_)[i] = x[i];
    }
    
    solver_->Iterate(options_.max_iterations,
                     options_.tolerance);
    
    lhs_->ExtractCopy(&x[0]);
}

void Aztec_Inverse_Operator::
check_class_invariants() const
{
    Assert(vector_operator_);
    Assert(vector_operator_->square());
}

string Aztec_Inverse_Operator::
description() const
{
    return "(Aztec_Inverse_Operator -> " + vector_operator_->description() + ")";
}
