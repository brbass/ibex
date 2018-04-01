#include "Meshless_Sweep.hh"

#include <iostream>
#if defined(ENABLE_OPENMP)
    #include <omp.h>
#else
    inline int omp_get_num_threads() {return 1;}
    inline int omp_get_thread_num() {return 0;}
#endif

#include "Amesos.h"
#include "AztecOO.h"
#include "AztecOO_ConditionNumber.h"
#include "BelosSolverFactory.hpp"
#include "BelosEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Vector.h"
#include "Ifpack.h"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_RCPStdSharedPtrConversions.hpp"
#include "Teuchos_ArrayRCP.hpp"

#include "Angular_Discretization.hh"
#include "Basis_Function.hh"
#include "Boundary_Source.hh"
#include "Cartesian_Plane.hh"
#include "Check.hh"
#include "Conversion.hh"
#include "Cross_Section.hh"
#include "Dimensional_Moments.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "Transport_Discretization.hh"
#include "XML_Node.hh"

using std::make_shared;
using std::pair;
using std::shared_ptr;
using std::string;
using std::vector;

Meshless_Sweep::
Meshless_Sweep(Options options,
               shared_ptr<Weak_Spatial_Discretization> spatial_discretization,
               shared_ptr<Angular_Discretization> angular_discretization,
               shared_ptr<Energy_Discretization> energy_discretization,
               shared_ptr<Transport_Discretization> transport_discretization):
    Sweep_Operator(Sweep_Type::ORDINATE,
                   transport_discretization),
    options_(options),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization)
{
}

void Meshless_Sweep::
initialize_solver()
{
    switch (options_.solver)
    {
    case Options::Solver::AMESOS:
        solver_ = make_shared<Amesos_Solver>(*this);
        break;
    case Options::Solver::AMESOS_PARALLEL:
        solver_ = make_shared<Amesos_Parallel_Solver>(*this);
        break;
    case Options::Solver::AZTEC:
        solver_ = make_shared<Aztec_Solver>(*this);
        break;
    case Options::Solver::AZTEC_IFPACK:
        solver_ = make_shared<Aztec_Ifpack_Solver>(*this);
        break;
    case Options::Solver::BELOS:
        solver_ = make_shared<Belos_Solver>(*this);
        break;
    case Options::Solver::BELOS_IFPACK:
        solver_ = make_shared<Belos_Ifpack_Solver>(*this);
        break;
    case Options::Solver::BELOS_IFPACK_RIGHT:
        solver_ = make_shared<Belos_Ifpack_Right_Solver>(*this);
        break;
    case Options::Solver::BELOS_IFPACK_RIGHT2:
        solver_ = make_shared<Belos_Ifpack_Right2_Solver>(*this);
        break;
    }
}

void Meshless_Sweep::
apply(vector<double> &x) const
{
    solver_->solve(x);
    update_augments(x);
}

void Meshless_Sweep::
update_augments(vector<double> &x) const
{
    int number_of_boundary_points = spatial_discretization_->number_of_boundary_points();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_groups = energy_discretization_->number_of_groups();
    int psi_size = transport_discretization_->psi_size();
    bool has_reflection = transport_discretization_->has_reflection();

    // No augments if there is no reflection
    if (!has_reflection)
    {
        return;
    }

    // Iterate over boundary indices
    for (int b = 0; b < number_of_boundary_points; ++b)
    {
        // Get global index
        int i = spatial_discretization_->boundary_basis(b)->index();

        for (int o = 0; o < number_of_ordinates; ++o)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                // Set boundary augment to current value of angular flux
                int k_b = psi_size + g + number_of_groups * (o + number_of_ordinates * b);
                int k_psi = g + number_of_groups * (o + number_of_ordinates * i);
                x[k_b] = x[k_psi];
            }
        }
    }
}

void Meshless_Sweep::
output(XML_Node output_node) const
{
    output_node.set_attribute(options_.solver_conversion()->convert(options_.solver),
                              "solver");
}

void Meshless_Sweep::
save_matrix_as_xml(int o,
                   int g,
                   XML_Node output_node) const
{
    // Create matrix node
    int number_of_points = spatial_discretization_->number_of_points();
    XML_Node matrix_node = output_node.append_child("matrix");
    matrix_node.set_attribute(o, "o");
    matrix_node.set_attribute(g, "g");
    matrix_node.set_attribute(number_of_points, "number_of_points");
    matrix_node.set_child_vector(spatial_discretization_->number_of_basis_functions(), "number_of_entries");
    
    for (int i = 0; i < number_of_points; ++i)
    {
        // Get row of matrix
        vector<int> indices;
        vector<double> values;
        get_matrix_row(i,
                       o,
                       g,
                       indices,
                       values);

        // Store row of matrix
        XML_Node row_node = matrix_node.append_child("row");
        row_node.set_attribute(i, "row_index");
        row_node.set_child_vector(indices,
                                  "column_indices");
        row_node.set_child_vector(values,
                                  "values");
    }
}

Meshless_Sweep::Sweep_Solver::
Sweep_Solver(Meshless_Sweep const &wrs):
    wrs_(wrs)
{
}

Meshless_Sweep::Trilinos_Solver::
Trilinos_Solver(Meshless_Sweep const &wrs):
    Sweep_Solver(wrs)
{
}

shared_ptr<Epetra_CrsMatrix> Meshless_Sweep::Trilinos_Solver::
get_matrix(int o,
           int g,
           shared_ptr<Epetra_Map> map) const
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    vector<int> const number_of_basis_functions = wrs_.spatial_discretization_->number_of_basis_functions();
    
    shared_ptr<Epetra_CrsMatrix> mat
        = make_shared<Epetra_CrsMatrix>(Copy, // Data access
                                        *map,
                                        &number_of_basis_functions[0], // Num entries per row
                                        true); // Static profile
    for (int i = 0; i < number_of_points; ++i)
    {
        vector<int> indices;
        vector<double> values;
        wrs_.get_matrix_row(i,
                            o,
                            g,
                            indices,
                            values);
        mat->InsertGlobalValues(i, // Row
                                number_of_basis_functions[i], // Num entries
                                &values[0],
                                &indices[0]);
    }
    mat->FillComplete();
    mat->OptimizeStorage();
    
    return mat;
}

shared_ptr<Epetra_CrsMatrix> Meshless_Sweep::Trilinos_Solver::
get_prec_matrix(shared_ptr<Epetra_Map> map) const
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    vector<int> const number_of_basis_functions = wrs_.spatial_discretization_->number_of_basis_functions();
    
    shared_ptr<Epetra_CrsMatrix> mat
        = make_shared<Epetra_CrsMatrix>(Copy, // Data access
                                        *map,
                                        &number_of_basis_functions[0], // Num entries per row
                                        true); // Static profile
    for (int i = 0; i < number_of_points; ++i)
    {
        vector<int> indices;
        vector<double> values;
        wrs_.get_prec_matrix_row(i,
                                 indices,
                                 values);
        mat->InsertGlobalValues(i, // Row
                                number_of_basis_functions[i], // Num entries
                                &values[0],
                                &indices[0]);
    }
    mat->FillComplete();
    mat->OptimizeStorage();
    
    return mat;
}

void Meshless_Sweep::Trilinos_Solver::
set_rhs(int o,
        int g,
        std::shared_ptr<Epetra_Vector> &rhs,
        vector<double> const &x) const
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    int number_of_groups = wrs_.energy_discretization_->number_of_groups();
    int k = g + number_of_groups * o;
    for (int i = 0; i < number_of_points; ++i)
    {
        double value;
        wrs_.get_rhs(i,
                     o,
                     g,
                     x,
                     value);
        
        (*rhs)[i] = value;
    }
}

void Meshless_Sweep::Trilinos_Solver::
check_aztec_convergence(shared_ptr<AztecOO> const solver) const
{
    bool converged;
    string message;
    double const *status = solver->GetAztecStatus();
    switch ((int) status[AZ_why])
    {
    case AZ_normal:
        converged = true;
        break;
    case AZ_param:
        converged = false;
        message = "AztecOO: Parameter not available";
        break;
    case AZ_breakdown:
        converged = false;
        message = "AztecOO: Numerical breakdown";
        break;
    case AZ_loss:
        converged = false;
        message = "AztecOO: Numerical loss of precision";
        break;
    case AZ_ill_cond:
        converged = false;
        message = "AztecOO: Ill-conditioned matrix";
        break;
    case AZ_maxits:
        converged = false;
        message = "AztecOO: Max iterations reached without convergence";
        break;
    }
    if (!converged)
    {
        if (wrs_.options_.quit_if_diverged)
        {
            AssertMsg(false, message);
        }
        else
        {
            std::cerr << message << std::endl;
        }
    }
}

Meshless_Sweep::Amesos_Solver::
Amesos_Solver(Meshless_Sweep const &wrs):
    Trilinos_Solver(wrs)
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    int number_of_groups = wrs_.energy_discretization_->number_of_groups();
    int number_of_ordinates = wrs_.angular_discretization_->number_of_ordinates();
    
    // Initialize communication
    comm_ = make_shared<Epetra_SerialComm>();
    map_ = make_shared<Epetra_Map>(number_of_points, 0, *comm_);
    
    // Initialize matrices and vectors
    lhs_ = make_shared<Epetra_Vector>(*map_);
    rhs_ = make_shared<Epetra_Vector>(*map_);
    lhs_->PutScalar(1.0);
    rhs_->PutScalar(1.0);
    mat_.resize(number_of_groups * number_of_ordinates);
    problem_.resize(number_of_groups * number_of_ordinates);
    solver_.resize(number_of_groups * number_of_ordinates);
    Amesos factory;
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k = g + number_of_groups * o;

            mat_[k] = get_matrix(o,
                                 g,
                                 map_);
            problem_[k]
                = make_shared<Epetra_LinearProblem>(mat_[k].get(),
                                                    lhs_.get(),
                                                    rhs_.get());
            
            solver_[k]
                = shared_ptr<Amesos_BaseSolver>(factory.Create("Klu",
                                                               *problem_[k]));
                 
            AssertMsg(solver_[k]->SymbolicFactorization() == 0, "Amesos solver symbolic factorization failed");
            AssertMsg(solver_[k]->NumericFactorization() == 0, "Amesos solver numeric factorization failed");
        }
    }
}

void Meshless_Sweep::Amesos_Solver::
solve(vector<double> &x) const
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    int number_of_groups = wrs_.energy_discretization_->number_of_groups();
    int number_of_ordinates = wrs_.angular_discretization_->number_of_ordinates();

    // Solve independently for each ordinate and group
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k = g + number_of_groups * o;
                
            // Set current RHS value
            set_rhs(o,
                    g,
                    rhs_,
                    x);
                
            // Solve, putting result into LHS
            AssertMsg(solver_[k]->Solve() == 0, "Amesos solver failed to solve");
            
            // Update solution value (overwrite x for this o and g)
            for (int i = 0; i < number_of_points; ++i)
            {
                int k_x = g + number_of_groups * (o + number_of_ordinates * i);
                x[k_x] = (*lhs_)[i];
            }
        }
    }
}

Meshless_Sweep::Amesos_Parallel_Solver::
Amesos_Parallel_Solver(Meshless_Sweep const &wrs):
    Trilinos_Solver(wrs)
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    int number_of_groups = wrs_.energy_discretization_->number_of_groups();
    int number_of_ordinates = wrs_.angular_discretization_->number_of_ordinates();
    
    // Initialize matrices
    comm_.resize(number_of_groups * number_of_ordinates);
    map_.resize(number_of_groups * number_of_ordinates);
    mat_.resize(number_of_groups * number_of_ordinates);
    lhs_.resize(number_of_groups * number_of_ordinates);
    rhs_.resize(number_of_groups * number_of_ordinates);
    problem_.resize(number_of_groups * number_of_ordinates);
    solver_.resize(number_of_groups * number_of_ordinates);
    Amesos factory;
    #pragma omp parallel for
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k = g + number_of_groups * o;

            comm_[k] = make_shared<Epetra_SerialComm>();
            map_[k] = make_shared<Epetra_Map>(number_of_points, 0, *comm_[k]);
            
            lhs_[k] = make_shared<Epetra_Vector>(*map_[k]);
            rhs_[k] = make_shared<Epetra_Vector>(*map_[k]);
            lhs_[k]->PutScalar(1.0);
            rhs_[k]->PutScalar(1.0);
            mat_[k] = get_matrix(o,
                                 g,
                                 map_[k]);
            problem_[k]
                = make_shared<Epetra_LinearProblem>(mat_[k].get(),
                                                    lhs_[k].get(),
                                                    rhs_[k].get());
            solver_[k]
                = shared_ptr<Amesos_BaseSolver>(factory.Create("Klu",
                                                               *problem_[k]));
                 
            AssertMsg(solver_[k]->SymbolicFactorization() == 0, "Amesos solver symbolic factorization failed");
            AssertMsg(solver_[k]->NumericFactorization() == 0, "Amesos solver numeric factorization failed");
        }
    }
}

void Meshless_Sweep::Amesos_Parallel_Solver::
solve(vector<double> &x) const
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    int number_of_groups = wrs_.energy_discretization_->number_of_groups();
    int number_of_ordinates = wrs_.angular_discretization_->number_of_ordinates();

    // Solve independently for each ordinate and group
    #pragma omp parallel for
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k = g + number_of_groups * o;
                
            // Set current RHS value
            set_rhs(o,
                    g,
                    rhs_[k],
                    x);
                
            // Solve, putting result into LHS
            AssertMsg(solver_[k]->Solve() == 0, "Amesos solver failed to solve");
            
            // Update solution value (overwrite x for this o and g)
            for (int i = 0; i < number_of_points; ++i)
            {
                int k_x = g + number_of_groups * (o + number_of_ordinates * i);
                x[k_x] = (*lhs_[k])[i];
            }
        }
    }
}

Meshless_Sweep::Aztec_Solver::
Aztec_Solver(Meshless_Sweep const &wrs):
    Trilinos_Solver(wrs)
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    comm_ = make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    map_ = make_shared<Epetra_Map>(number_of_points, 0, *comm_);
    lhs_ = make_shared<Epetra_Vector>(*map_);
    rhs_ = make_shared<Epetra_Vector>(*map_);
    lhs_->PutScalar(1.0);
    rhs_->PutScalar(1.0);
}

void Meshless_Sweep::Aztec_Solver::
solve(vector<double> &x) const
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    int number_of_groups = wrs_.energy_discretization_->number_of_groups();
    int number_of_ordinates = wrs_.angular_discretization_->number_of_ordinates();

    // Solve independently for each ordinate and group
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            // Set current RHS value
            set_rhs(o,
                    g,
                    rhs_,
                    x);

            // Get matrix
            std::shared_ptr<Epetra_CrsMatrix> mat = get_matrix(o,
                                                               g,
                                                               map_);

            // Get linear problem
            int k = g + number_of_groups * o;
            std::shared_ptr<Epetra_LinearProblem> problem
                = make_shared<Epetra_LinearProblem>(mat.get(),
                                                    lhs_.get(),
                                                    rhs_.get());
            
            // Get solver
            shared_ptr<AztecOO> solver
                = make_shared<AztecOO>(*problem);
            solver->SetAztecOption(AZ_solver, AZ_gmres);
            solver->SetAztecOption(AZ_kspace, wrs_.options_.kspace);
            solver->SetAztecOption(AZ_precond, AZ_none);
            if (wrs_.options_.print)
            {
                solver->SetAztecOption(AZ_output, AZ_all);
            }
            else
            {
                solver->SetAztecOption(AZ_output, AZ_warnings);
            }
            
            // Solve, putting result into LHS
            solver->Iterate(wrs_.options_.max_iterations,
                            wrs_.options_.tolerance);

            // Check to ensure solver converged
            check_aztec_convergence(solver);
            
            // Update solution value (overwrite x for this o and g)
            for (int i = 0; i < number_of_points; ++i)
            {
                int k_x = g + number_of_groups * (o + number_of_ordinates * i);
                x[k_x] = (*lhs_)[i];
            }
        }
    }
}

Meshless_Sweep::Aztec_Ifpack_Solver::
Aztec_Ifpack_Solver(Meshless_Sweep const &wrs):
    Trilinos_Solver(wrs)
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    int number_of_groups = wrs_.energy_discretization_->number_of_groups();
    int number_of_ordinates = wrs_.angular_discretization_->number_of_ordinates();
    
    // Initialize communication
    comm_ = make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    map_ = make_shared<Epetra_Map>(number_of_points, 0, *comm_);
    
    // Initialize matrices
    lhs_ = make_shared<Epetra_Vector>(*map_);
    rhs_ = make_shared<Epetra_Vector>(*map_);
    lhs_->PutScalar(1.0);
    rhs_->PutScalar(1.0);
    mat_.resize(number_of_groups * number_of_ordinates);
    problem_.resize(number_of_groups * number_of_ordinates);
    if (wrs_.options_.use_preconditioner)
    {
        prec_.resize(number_of_groups * number_of_ordinates);
    }
    solver_.resize(number_of_groups * number_of_ordinates);
    
    Ifpack ifp_factory;
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k = g + number_of_groups * o;

            // Get matrix and problem
            mat_[k] = get_matrix(o,
                                 g,
                                 map_);
            problem_[k]
                = make_shared<Epetra_LinearProblem>(mat_[k].get(),
                                                    lhs_.get(),
                                                    rhs_.get());

            if (wrs_.options_.use_preconditioner)
            {
                // Create preconditioner
                // ILU requires an int "fact: level-of-fill"
                // ILUT requires a double "fact: ilut level-of-fill"
                prec_[k]
                    = shared_ptr<Ifpack_Preconditioner>(ifp_factory.Create("ILUT",
                                                                           mat_[k].get()));
                Teuchos::ParameterList prec_list;
                prec_list.set("fact: drop tolerance", wrs_.options_.drop_tolerance);
                // prec_list.set("fact: level-of-fill", wrs_.options_.level_of_fill);
                prec_list.set("fact: ilut level-of-fill", wrs_.options_.level_of_fill);
                prec_[k]->SetParameters(prec_list);
                prec_[k]->Initialize();
                prec_[k]->Compute();
                Assert(prec_[k]->IsInitialized() == true);
                Assert(prec_[k]->IsComputed() == true);
            }
            
            // Initialize solver
            solver_[k] = make_shared<AztecOO>(*problem_[k]);
            solver_[k]->SetAztecOption(AZ_solver, AZ_gmres);
            solver_[k]->SetAztecOption(AZ_kspace, wrs_.options_.kspace);
            if (wrs_.options_.use_preconditioner)
            {
                solver_[k]->SetPrecOperator(prec_[k].get());
            }
            if (wrs_.options_.print)
            {
                solver_[k]->SetAztecOption(AZ_output, AZ_all);
            }
            else
            {
                solver_[k]->SetAztecOption(AZ_output, AZ_warnings);
            }
        }
    }
}

void Meshless_Sweep::Aztec_Ifpack_Solver::
solve(vector<double> &x) const
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    int number_of_groups = wrs_.energy_discretization_->number_of_groups();
    int number_of_ordinates = wrs_.angular_discretization_->number_of_ordinates();

    // Solve independently for each ordinate and group
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k = g + number_of_groups * o;

            // Set the LHS value
            // If not set, the problem can converge too quickly
            // and the solver thinks it's incorrectly converged
            lhs_->PutScalar(1.0);
            
            // Set current RHS value
            set_rhs(o,
                    g,
                    rhs_,
                    x);
            
            // Solve, putting result into LHS
            solver_[k]->Iterate(wrs_.options_.max_iterations,
                                wrs_.options_.tolerance);
            
            // Check to ensure solver converged
            check_aztec_convergence(solver_[k]);
            
            // Update solution value (overwrite x for this o and g)
            for (int i = 0; i < number_of_points; ++i)
            {
                int k_x = g + number_of_groups * (o + number_of_ordinates * i);
                x[k_x] = (*lhs_)[i];
            }
        }
    }
}

Meshless_Sweep::Belos_Solver::
Belos_Solver(Meshless_Sweep const &wrs):
    Trilinos_Solver(wrs)
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    
    #pragma omp parallel
    {
        int number_of_threads = omp_get_num_threads();
        int t = omp_get_thread_num();

        #pragma omp single
        {
            // Initialize data pointers
            comm_.resize(number_of_threads);
            map_.resize(number_of_threads);
            lhs_.resize(number_of_threads);
            rhs_.resize(number_of_threads);
            problem_.resize(number_of_threads);
            solver_.resize(number_of_threads);
        }

        // Get comm and map
        comm_[t] = make_shared<Epetra_SerialComm>();
        map_[t] = make_shared<Epetra_Map>(number_of_points, 0, *comm_[t]);
        
        // Get vectors
        lhs_[t] = make_shared<Epetra_Vector>(*map_[t]);
        rhs_[t] = make_shared<Epetra_Vector>(*map_[t]);
        lhs_[t]->PutScalar(1.0);
        rhs_[t]->PutScalar(1.0);

        // Get problem and solver
        shared_ptr<Teuchos::ParameterList> belos_list
            = make_shared<Teuchos::ParameterList>();
        belos_list->set("Num Blocks", wrs_.options_.kspace);
        belos_list->set("Maximum Iterations", wrs_.options_.max_iterations);
        belos_list->set("Maximum Restarts", wrs_.options_.max_restarts);
        belos_list->set("Convergence Tolerance", wrs_.options_.tolerance);
        if (wrs_.options_.print)
        {
            belos_list->set("Verbosity", Belos::IterationDetails + Belos::TimingDetails + Belos::FinalSummary);
        }
        else
        {
            belos_list->set("Verbosity", Belos::Errors + Belos::Warnings);
        }
        #pragma omp critical
        {
            problem_[t] = make_shared<BelosLinearProblem>();
            problem_[t]->setLHS(Teuchos::rcp(lhs_[t]));
            problem_[t]->setRHS(Teuchos::rcp(rhs_[t]));
            solver_[t] = make_shared<BelosSolver>(Teuchos::rcp(problem_[t]),
                                                  Teuchos::rcp(belos_list));
        }
    }
}

void Meshless_Sweep::Belos_Solver::
solve(vector<double> &x) const
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    int number_of_groups = wrs_.energy_discretization_->number_of_groups();
    int number_of_ordinates = wrs_.angular_discretization_->number_of_ordinates();

    // Solve independently for each ordinate and group
    #pragma omp parallel
    {
        int number_of_threads = omp_get_num_threads();
        int t = omp_get_thread_num();
        Assert(problem_.size() == number_of_threads);
        
        #pragma omp for
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k = g + number_of_groups * o;
                string description = std::to_string(o) + "_" + std::to_string(g);
                
                // Set current RHS value
                set_rhs(o,
                        g,
                        rhs_[t],
                        x);
                
                // Initialize LHS to 1.0 to avoid implicit residual problems
                lhs_[t]->PutScalar(1.0);

                // Get matrix
                shared_ptr<Epetra_CrsMatrix> mat
                    = get_matrix(o,
                                 g,
                                 map_[t]);

                // Set up problem
                problem_[t]->setOperator(Teuchos::rcp(mat));
                AssertMsg(problem_[t]->setProblem(), description);
            
                // Solve, putting result into LHS
                try
                {
                    Belos::ReturnType belos_result
                        = solver_[t]->solve();
                
                    if (wrs_.options_.quit_if_diverged)
                    {
                        AssertMsg(belos_result == Belos::Converged, description);
                    }
                }
                catch (Belos::StatusTestError const &error)
                {
                    AssertMsg(false, "Belos status test failed, " + description);
                }
                // std::cout << solver_[k]->getNumIters() << std::endl;
            
                // Update solution value (overwrite x for this o and g)
                for (int i = 0; i < number_of_points; ++i)
                {
                    int k_x = g + number_of_groups * (o + number_of_ordinates * i);
                    x[k_x] = (*lhs_[t])[i];
                }
            }
        }

    }
}

Meshless_Sweep::Belos_Ifpack_Solver::
Belos_Ifpack_Solver(Meshless_Sweep const &wrs):
    Trilinos_Solver(wrs)
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    int number_of_groups = wrs_.energy_discretization_->number_of_groups();
    int number_of_ordinates = wrs_.angular_discretization_->number_of_ordinates();
    
    // Initialize matrices
    comm_.resize(number_of_groups * number_of_ordinates);
    map_.resize(number_of_groups * number_of_ordinates);
    lhs_.resize(number_of_groups * number_of_ordinates);
    rhs_.resize(number_of_groups * number_of_ordinates);
    mat_.resize(number_of_groups * number_of_ordinates);
    if (wrs_.options_.use_preconditioner)
    {
        prec_.resize(number_of_groups * number_of_ordinates);
    }
    problem_.resize(number_of_groups * number_of_ordinates);
    solver_.resize(number_of_groups * number_of_ordinates);
    
    #pragma omp parallel for
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k = g + number_of_groups * o;
            string description = std::to_string(o) + "_" + std::to_string(g);

            // Get comm and map
            comm_[k] = make_shared<Epetra_SerialComm>();
            map_[k] = make_shared<Epetra_Map>(number_of_points, 0, *comm_[k]);

            // Get vectors and matrix
            lhs_[k] = make_shared<Epetra_Vector>(*map_[k]);
            rhs_[k] = make_shared<Epetra_Vector>(*map_[k]);
            lhs_[k]->PutScalar(1.0);
            rhs_[k]->PutScalar(1.0);
            mat_[k] = get_matrix(o,
                                 g,
                                 map_[k]);

            // Get preconditioner
            if (wrs_.options_.use_preconditioner)
            {
                Ifpack factory;
                shared_ptr<Ifpack_Preconditioner> temp_prec
                    = shared_ptr<Ifpack_Preconditioner>(factory.Create("ILUT",
                                                                       mat_[k].get()));
                Teuchos::ParameterList prec_list;
                prec_list.set("fact: drop tolerance", wrs_.options_.drop_tolerance);
                prec_list.set("fact: ilut level-of-fill", wrs_.options_.level_of_fill);
                temp_prec->SetParameters(prec_list);
                temp_prec->Initialize();
                temp_prec->Compute();
                AssertMsg(temp_prec->IsInitialized() == true, description);
                AssertMsg(temp_prec->IsComputed() == true, description);
                
                prec_[k]
                    = make_shared<BelosPreconditioner>(Teuchos::rcp(temp_prec));
            }

            // Get problem
            #pragma omp critical
            {
                problem_[k]
                    = make_shared<BelosLinearProblem>(Teuchos::rcp(mat_[k]),
                                                      Teuchos::rcp(lhs_[k]),
                                                      Teuchos::rcp(rhs_[k]));
                if (wrs_.options_.use_preconditioner)
                {
                    problem_[k]->setLeftPrec(Teuchos::rcp(prec_[k]));
                }
                AssertMsg(problem_[k]->setProblem(), description);
            }
            
            // Get solver
            shared_ptr<Teuchos::ParameterList> belos_list
                = make_shared<Teuchos::ParameterList>();
            belos_list->set("Num Blocks", wrs_.options_.kspace);
            belos_list->set("Maximum Iterations", wrs_.options_.max_iterations);
            belos_list->set("Maximum Restarts", wrs_.options_.max_restarts);
            belos_list->set("Convergence Tolerance", wrs_.options_.tolerance);
            if (wrs_.options_.print)
            {
                belos_list->set("Verbosity", Belos::IterationDetails + Belos::TimingDetails + Belos::FinalSummary);
            }
            else
            {
                belos_list->set("Verbosity", Belos::Errors + Belos::Warnings);
            }
            #pragma omp critical
            {
                solver_[k]
                    = make_shared<BelosSolver>(Teuchos::rcp(problem_[k]),
                                               Teuchos::rcp(belos_list));
            }
        }
    }
}

void Meshless_Sweep::Belos_Ifpack_Solver::
solve(vector<double> &x) const
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    int number_of_groups = wrs_.energy_discretization_->number_of_groups();
    int number_of_ordinates = wrs_.angular_discretization_->number_of_ordinates();

    // Solve independently for each ordinate and group
    #pragma omp parallel for
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k = g + number_of_groups * o;
            string description = std::to_string(o) + "_" + std::to_string(g);
                
            // Set current RHS value
            set_rhs(o,
                    g,
                    rhs_[k],
                    x);

            // Initialize LHS to 1.0 to avoid implicit residual problems
            lhs_[k]->PutScalar(1.0);

            // Set up problem
            AssertMsg(problem_[k]->setProblem(), description);
            
            // Solve, putting result into LHS
            try
            {
                Belos::ReturnType belos_result
                    = solver_[k]->solve();
                
                if (wrs_.options_.quit_if_diverged)
                {
                    AssertMsg(belos_result == Belos::Converged, description);
                }
            }
            catch (Belos::StatusTestError const &error)
            {
                AssertMsg(false, "Belos status test failed, " + description);
            }
            // std::cout << solver_[k]->getNumIters() << std::endl;
            
            // Update solution value (overwrite x for this o and g)
            for (int i = 0; i < number_of_points; ++i)
            {
                int k_x = g + number_of_groups * (o + number_of_ordinates * i);
                x[k_x] = (*lhs_[k])[i];
            }
        }
    }
}

Meshless_Sweep::Belos_Ifpack_Right_Solver::
Belos_Ifpack_Right_Solver(Meshless_Sweep const &wrs):
    Trilinos_Solver(wrs)
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    
    #pragma omp parallel
    {
        int number_of_threads = omp_get_num_threads();
        int t = omp_get_thread_num();

        #pragma omp single
        {
            // Initialize data pointers
            comm_.resize(number_of_threads);
            map_.resize(number_of_threads);
            lhs_.resize(number_of_threads);
            rhs_.resize(number_of_threads);
            prec_.resize(number_of_threads);
            prec_mat_.resize(number_of_threads);
            problem_.resize(number_of_threads);
            solver_.resize(number_of_threads);
        }
        
        // Get comm and map
        comm_[t] = make_shared<Epetra_SerialComm>();
        map_[t] = make_shared<Epetra_Map>(number_of_points, 0, *comm_[t]);
        
        // Get vectors
        lhs_[t] = make_shared<Epetra_Vector>(*map_[t]);
        rhs_[t] = make_shared<Epetra_Vector>(*map_[t]);
        lhs_[t]->PutScalar(1.0);
        rhs_[t]->PutScalar(1.0);
        
        // Get preconditioner
        if (wrs_.options_.use_preconditioner)
        {
            prec_mat_[t] = get_prec_matrix(map_[t]);

            Ifpack factory;
            shared_ptr<Ifpack_Preconditioner> temp_prec
                = shared_ptr<Ifpack_Preconditioner>(factory.Create("ILUT",
                                                                   prec_mat_[t].get()));
            Teuchos::ParameterList prec_list;
            prec_list.set("fact: drop tolerance", wrs_.options_.drop_tolerance);
            prec_list.set("fact: ilut level-of-fill", wrs_.options_.level_of_fill);
            temp_prec->SetParameters(prec_list);
            temp_prec->Initialize();
            temp_prec->Compute();
            AssertMsg(temp_prec->IsInitialized() == true, std::to_string(t));
            AssertMsg(temp_prec->IsComputed() == true, std::to_string(t));
                
            prec_[t]
                = make_shared<BelosPreconditioner>(Teuchos::rcp(temp_prec));
        }
        
        // Get problem and solver
        shared_ptr<Teuchos::ParameterList> belos_list
            = make_shared<Teuchos::ParameterList>();
        belos_list->set("Num Blocks", wrs_.options_.kspace);
        belos_list->set("Maximum Iterations", wrs_.options_.max_iterations);
        belos_list->set("Maximum Restarts", wrs_.options_.max_restarts);
        belos_list->set("Convergence Tolerance", wrs_.options_.tolerance);
        if (wrs_.options_.print)
        {
            belos_list->set("Verbosity", Belos::IterationDetails + Belos::TimingDetails + Belos::FinalSummary);
        }
        else
        {
            belos_list->set("Verbosity", Belos::Errors + Belos::Warnings);
        }
        #pragma omp critical
        {
            problem_[t] = make_shared<BelosLinearProblem>();
            if (wrs_.options_.use_preconditioner)
            {
                if (wrs_.options_.force_left)
                {
                    problem_[t]->setLeftPrec(Teuchos::rcp(prec_[t]));
                }
                else
                {
                    problem_[t]->setRightPrec(Teuchos::rcp(prec_[t]));
                }
            }
            problem_[t]->setLHS(Teuchos::rcp(lhs_[t]));
            problem_[t]->setRHS(Teuchos::rcp(rhs_[t]));
            solver_[t] = make_shared<BelosSolver>(Teuchos::rcp(problem_[t]),
                                                  Teuchos::rcp(belos_list));
        }
    }
}

void Meshless_Sweep::Belos_Ifpack_Right_Solver::
solve(vector<double> &x) const
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    int number_of_groups = wrs_.energy_discretization_->number_of_groups();
    int number_of_ordinates = wrs_.angular_discretization_->number_of_ordinates();

    // Solve independently for each ordinate and group
    #pragma omp parallel
    {
        int number_of_threads = omp_get_num_threads();
        int t = omp_get_thread_num();
        Assert(problem_.size() == number_of_threads);
        
        #pragma omp for
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k = g + number_of_groups * o;
                string description = std::to_string(o) + "_" + std::to_string(g);
                
                // Set current RHS value
                set_rhs(o,
                        g,
                        rhs_[t],
                        x);
                
                // Initialize LHS to 1.0 to avoid implicit residual problems
                lhs_[t]->PutScalar(1.0);

                // Get matrix
                shared_ptr<Epetra_CrsMatrix> mat
                    = get_matrix(o,
                                 g,
                                 map_[t]);
                
                // Set up problem
                problem_[t]->setOperator(Teuchos::rcp(mat));
                AssertMsg(problem_[t]->setProblem(), description);
                
                // Solve, putting result into LHS
                try
                {
                    Belos::ReturnType belos_result
                        = solver_[t]->solve();
                
                    if (wrs_.options_.quit_if_diverged)
                    {
                        AssertMsg(belos_result == Belos::Converged, description);
                    }
                }
                catch (Belos::StatusTestError const &error)
                {
                    AssertMsg(false, "Belos status test failed, " + description);
                }
                // std::cout << solver_[k]->getNumIters() << std::endl;
            
                // Update solution value (overwrite x for this o and g)
                for (int i = 0; i < number_of_points; ++i)
                {
                    int k_x = g + number_of_groups * (o + number_of_ordinates * i);
                    x[k_x] = (*lhs_[t])[i];
                }
            }
        }

    }
}

Meshless_Sweep::Belos_Ifpack_Right2_Solver::
Belos_Ifpack_Right2_Solver(Meshless_Sweep const &wrs):
    Trilinos_Solver(wrs)
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    int number_of_groups = wrs_.energy_discretization_->number_of_groups();
    int number_of_ordinates = wrs_.angular_discretization_->number_of_ordinates();

    #pragma omp parallel
    {
        int number_of_threads = omp_get_num_threads();
        int t = omp_get_thread_num();

        #pragma omp single
        {
            // Initialize data pointers
            comm_.resize(number_of_threads);
            map_.resize(number_of_threads);
            lhs_.resize(number_of_threads);
            rhs_.resize(number_of_threads);
            prec_.resize(number_of_threads);
            prec_mat_.resize(number_of_threads);
            mat_.resize(number_of_groups * number_of_ordinates);
            problem_.resize(number_of_groups * number_of_ordinates);
            solver_.resize(number_of_groups * number_of_ordinates);
        }
        
        // Get comm and map
        comm_[t] = make_shared<Epetra_SerialComm>();
        map_[t] = make_shared<Epetra_Map>(number_of_points, 0, *comm_[t]);
        
        // Get vectors
        lhs_[t] = make_shared<Epetra_Vector>(*map_[t]);
        rhs_[t] = make_shared<Epetra_Vector>(*map_[t]);
        lhs_[t]->PutScalar(1.0);
        rhs_[t]->PutScalar(1.0);
        
        // Get preconditioner
        if (wrs_.options_.use_preconditioner)
        {
            prec_mat_[t] = get_prec_matrix(map_[t]);

            Ifpack factory;
            shared_ptr<Ifpack_Preconditioner> temp_prec
                = shared_ptr<Ifpack_Preconditioner>(factory.Create("ILUT",
                                                                   prec_mat_[t].get()));
            Teuchos::ParameterList prec_list;
            prec_list.set("fact: drop tolerance", wrs_.options_.drop_tolerance);
            prec_list.set("fact: ilut level-of-fill", wrs_.options_.level_of_fill);
            temp_prec->SetParameters(prec_list);
            temp_prec->Initialize();
            temp_prec->Compute();
            AssertMsg(temp_prec->IsInitialized() == true, std::to_string(t));
            AssertMsg(temp_prec->IsComputed() == true, std::to_string(t));
                
            prec_[t]
                = make_shared<BelosPreconditioner>(Teuchos::rcp(temp_prec));
        }
        
        // Get problem and solver
        shared_ptr<Teuchos::ParameterList> belos_list
            = make_shared<Teuchos::ParameterList>();
        belos_list->set("Num Blocks", wrs_.options_.kspace);
        belos_list->set("Maximum Iterations", wrs_.options_.max_iterations);
        belos_list->set("Maximum Restarts", wrs_.options_.max_restarts);
        belos_list->set("Convergence Tolerance", wrs_.options_.tolerance);
        if (wrs_.options_.print)
        {
            belos_list->set("Verbosity", Belos::IterationDetails + Belos::TimingDetails + Belos::FinalSummary);
        }
        else
        {
            belos_list->set("Verbosity", Belos::Errors + Belos::Warnings);
        }
        
        #pragma omp for schedule(static)
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                string description = std::to_string(o) + "_" + std::to_string(g);
                int k = g + number_of_groups * o;
                
                mat_[k] = get_matrix(o,
                                     g,
                                     map_[t]);
                
                #pragma omp critical
                {
                    problem_[k] = make_shared<BelosLinearProblem>(Teuchos::rcp(mat_[k]),
                                                                  Teuchos::rcp(lhs_[t]),
                                                                  Teuchos::rcp(rhs_[t]));
                    if (wrs_.options_.use_preconditioner)
                    {
                        if (wrs_.options_.force_left)
                        {
                            problem_[k]->setLeftPrec(Teuchos::rcp(prec_[t]));
                        }
                        else
                        {
                            problem_[k]->setRightPrec(Teuchos::rcp(prec_[t]));
                        }
                    }
                    AssertMsg(problem_[k]->setProblem(), description);
                    
                    solver_[k] = make_shared<BelosSolver>(Teuchos::rcp(problem_[k]),
                                                          Teuchos::rcp(belos_list));
                }
            }
        }
    }
}

void Meshless_Sweep::Belos_Ifpack_Right2_Solver::
solve(vector<double> &x) const
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    int number_of_groups = wrs_.energy_discretization_->number_of_groups();
    int number_of_ordinates = wrs_.angular_discretization_->number_of_ordinates();

    // Solve independently for each ordinate and group
    #pragma omp parallel
    {
        int number_of_threads = omp_get_num_threads();
        int t = omp_get_thread_num();
        
        #pragma omp for schedule(static)
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k = g + number_of_groups * o;
                string description = std::to_string(o) + "_" + std::to_string(g);
                
                // Set current RHS value
                set_rhs(o,
                        g,
                        rhs_[t],
                        x);
                
                // Initialize LHS to 1.0 to avoid implicit residual problems
                lhs_[t]->PutScalar(1.0);

                // Set up problem
                AssertMsg(problem_[k]->setProblem(), description);
                
                // Solve, putting result into LHS
                try
                {
                    Belos::ReturnType belos_result
                        = solver_[k]->solve();
                
                    if (wrs_.options_.quit_if_diverged)
                    {
                        AssertMsg(belos_result == Belos::Converged, description);
                    }
                }
                catch (Belos::StatusTestError const &error)
                {
                    AssertMsg(false, "Belos status test failed, " + description);
                }

                // std::cout << solver_[k]->getNumIters() << std::endl;
            
                // Update solution value (overwrite x for this o and g)
                for (int i = 0; i < number_of_points; ++i)
                {
                    int k_x = g + number_of_groups * (o + number_of_ordinates * i);
                    x[k_x] = (*lhs_[t])[i];
                }
            }
        }
    }
}

shared_ptr<Conversion<Meshless_Sweep::Options::Solver, string> > Meshless_Sweep::Options::
solver_conversion() const
{
    vector<pair<Solver, string> > conversions
        = {{Solver::AMESOS, "amesos"},
           {Solver::AMESOS_PARALLEL, "amesos_parallel"},
           {Solver::AZTEC, "aztec"},
           {Solver::AZTEC_IFPACK, "aztec_ifpack"},
           {Solver::BELOS, "belos"},
           {Solver::BELOS_IFPACK, "belos_ifpack"},
           {Solver::BELOS_IFPACK_RIGHT, "belos_ifpack_right"},
           {Solver::BELOS_IFPACK_RIGHT2, "belos_ifpack_right2"}};
           
    return make_shared<Conversion<Solver, string> >(conversions);
}
