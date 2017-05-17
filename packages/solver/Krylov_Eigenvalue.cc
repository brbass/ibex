#include "Krylov_Eigenvalue.hh"

#include <AnasaziBlockKrylovSchurSolMgr.hpp>
#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziEpetraAdapter.hpp>
#include <AnasaziGeneralizedDavidsonSolMgr.hpp>
#include <AztecOO.h>
#include <Epetra_MultiVector.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>
#include "Teuchos_RCPStdSharedPtrConversions.hpp"

#include "Angular_Discretization.hh"
#include "Aztec_Inverse_Operator.hh"
#include "Convergence_Measure.hh"
#include "Energy_Discretization.hh"
#include "Epetra_Operator_Interface.hh"
#include "Spatial_Discretization.hh"
#include "Transport_Discretization.hh"
#include "Vector_Operator.hh"
#include "Vector_Operator_Functions.hh"
#include "XML_Node.hh"

typedef Anasazi::BasicEigenproblem<double, Epetra_MultiVector, Epetra_Operator> Anasazi_Eigenproblem;
typedef Anasazi::GeneralizedDavidsonSolMgr<double, Epetra_MultiVector, Epetra_Operator> Anasazi_SolverManager;
typedef Anasazi::Eigensolution<double, Epetra_MultiVector> Anasazi_Eigensolution;

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ParameterList;
using Teuchos::parameterList;
using Teuchos::get_shared_ptr;

Krylov_Eigenvalue::
Krylov_Eigenvalue(Options options,
                    shared_ptr<Spatial_Discretization> spatial_discretization,
                    shared_ptr<Angular_Discretization> angular_discretization,
                    shared_ptr<Energy_Discretization> energy_discretization,
                    shared_ptr<Transport_Discretization> transport_discretization,
                  shared_ptr<Vector_Operator> fission_operator,
                  shared_ptr<Vector_Operator> flux_operator,
                  vector<shared_ptr<Vector_Operator> > value_operators):
    Solver(options.solver_print,
           Solver::Type::K_EIGENVALUE),
    options_(options),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    transport_discretization_(transport_discretization),
    fission_operator_(fission_operator),
    flux_operator_(flux_operator),
    value_operators_(value_operators)
{
    check_class_invariants();
}

void Krylov_Eigenvalue::
solve()
{
    int phi_size = transport_discretization_->phi_size();
    int number_of_augments = transport_discretization_->number_of_augments();

    // Initialize result
    result_ = make_shared<Result>();

    // Initialize comm and map
    shared_ptr<Epetra_Comm> comm = make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    shared_ptr<Epetra_Map> map = make_shared<Epetra_Map>(phi_size + number_of_augments, 0, *comm);
    
    // Get inverse operator
    Aztec_Inverse_Operator::Options options;
    options.max_iterations = options_.max_inverse_iterations;
    options.kspace = options_.kspace;
    options.solver_print = options_.solver_print;
    options.tolerance = options_.tolerance;
    shared_ptr<Inverse_Operator> inverse_operator
        = make_shared<Aztec_Inverse_Operator>(options,
                                              flux_operator_,
                                              comm,
                                              map);
    
    // Get combined operator
    shared_ptr<Vector_Operator> full_operator
        = inverse_operator * fission_operator_;
    
    // Get Epetra operator
    shared_ptr<Epetra_Operator> oper
        = make_shared<Epetra_Operator_Interface>(comm,
                                                 map,
                                                 full_operator);
    
    // Create other Trilinos members
    shared_ptr<Epetra_MultiVector> init
        = make_shared<Epetra_MultiVector>(*map,
                                          options_.block_size);
    init->PutScalar(1.0);

    // Create problem
    shared_ptr<Anasazi_Eigenproblem> problem
        = make_shared<Anasazi_Eigenproblem>();
    problem->setOperator(rcp(oper));
    problem->setInitVec(rcp(init));
    problem->setNEV(options_.number_of_eigenvalues);
    problem->setHermitian(false);
    bool initialization_successful = problem->setProblem();
    Assert(initialization_successful);

    // Set general eigenvalue parameters
    shared_ptr<ParameterList> params
        = get_shared_ptr(parameterList());
    params->set("Maximum Iterations", options_.max_iterations);
    params->set("Block Size", options_.block_size);
    params->set("Convergence Tolerance", options_.tolerance);
    params->set("Which", "LR");
    params->set("Maximum Restarts", options_.max_iterations);
    params->set("Relative Convergence Tolerance", true);
    if (options_.solver_print)
    {
        int verbosity = Anasazi::IterationDetails + Anasazi::TimingDetails + Anasazi::FinalSummary;
        params->set("Verbosity", verbosity);
        params->set("Output Frequency", 1);
    }
    else
    {
        int verbosity = Anasazi::Errors;
        params->set("Verbosity", verbosity);
    }
    
    // Set GeneralizedDavidson parameters
    int kspace = min(options_.kspace, phi_size + number_of_augments);
    params->set("Maximum Subspace Dimension", kspace);
    params->set("Restart Dimension", int(ceil(kspace / 3)));
    params->set("Initial Guess", "User");
    
    // Get solver
    shared_ptr<Anasazi_SolverManager> solver
        = make_shared<Anasazi_SolverManager>(rcp(problem), *params);

    // Solve problem
    Anasazi::ReturnType converged = solver->solve();
    result_->total_iterations = solver->getNumIters();

    // Check for convergence
    if (converged != Anasazi::ReturnType::Converged)
    {
        cerr << "Eigenvalue solve did not converge" << endl;
    }

    // Get solution
    shared_ptr<Anasazi_Eigensolution const> solution
        = make_shared<Anasazi_Eigensolution> (problem->getSolution());
    
    // Get eigenvalue
    if (solution->numVecs < 1)
    {
        cerr << "No eigenvalues computed" << endl;
        result_->k_eigenvalue = -1;
        return;
    }
    else
    {
        result_->k_eigenvalue = solution->Evals[0].realpart;
    }

    // Get eigenvector
    shared_ptr<Epetra_MultiVector> eigenvectors
        = get_shared_ptr(solution->Evecs);
    vector<double> &coefficients = result_->coefficients;
    coefficients.resize(eigenvectors->Stride() * eigenvectors->NumVectors());
    if (eigenvectors->NumVectors() > 1)
    {
        cout << "Discarding additional eigenvectors" << endl;
    }
    eigenvectors->ExtractCopy(&coefficients[0],
                              eigenvectors->Stride());
    coefficients.resize(phi_size);

    // Get flux
    int number_of_values = value_operators_.size();
    result_->phi.resize(number_of_values);
    for (int i = 0; i < number_of_values; ++i)
    {
        vector<double> &phi = result_->phi[i];
        phi = coefficients;
        (*value_operators_[i])(phi);
    }
}

void Krylov_Eigenvalue::
output(XML_Node output_node) const
{
    // Output options
    output_node.set_attribute(options_.max_inverse_iterations,
                              "max_inverse_iterations");
    output_node.set_attribute(options_.max_iterations,
                               "max_iterations");
    output_node.set_attribute(options_.kspace,
                              "kspace");
    output_node.set_attribute(options_.solver_print,
                              "solver_print");
    output_node.set_attribute(options_.tolerance,
                              "tolerance");
    
    // Output results
    output_result(output_node,
                  result_);
}

void Krylov_Eigenvalue::
check_class_invariants() const
{
    Assert(spatial_discretization_);
    Assert(angular_discretization_);
    Assert(energy_discretization_);
    Assert(transport_discretization_);
    Assert(fission_operator_);
    Assert(flux_operator_);
    for (std::shared_ptr<Vector_Operator> oper : value_operators_)
    {
        Assert(oper);
    }
}
