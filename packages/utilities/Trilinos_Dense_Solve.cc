#include "Trilinos_Dense_Solve.hh"

#include <string>
#include <vector>

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseSolver.h>

#include <Amesos.h>
#include <AztecOO.h>
#ifdef EPETRA_MPI
#  include <mpi.h>
#  include <Epetra_MpiComm.h>
#else
#  include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>

#include "Check.hh"

using namespace std;

Trilinos_Dense_Solve::
Trilinos_Dense_Solve(int size,
                     Solver solver):
    size_(size),
    solver_(solver)
{
}

void Trilinos_Dense_Solve::
solve(vector<double> &a_data,
      vector<double> &b_data,
      vector<double> &x_data) const
{
    switch (solver_)
    {
    case Solver::EPETRA:
        return epetra_solve(a_data,
                            b_data,
                            x_data);
    case Solver::AMESOS:
        return amesos_solve(a_data,
                            b_data,
                            x_data);
    case Solver::AZTEC:
        return aztec_solve(a_data,
                           b_data,
                           x_data);
    }
}

// Iterative matrix solution: faster than GSL for n>300
void Trilinos_Dense_Solve::
epetra_solve(vector<double> &a_data,
             vector<double> &b_data,
             vector<double> &x_data) const
{
    Check(a_data.size() == size_ * size_);
    Check(b_data.size() == size_);
    Check(x_data.size() == size_);

    Epetra_SerialDenseMatrix a(View, &a_data[0], size_, size_, size_);
    Epetra_SerialDenseVector b(View, &b_data[0], size_);
    Epetra_SerialDenseVector x(size_);

    // cout << a << endl;
    // cout << b << endl;
    // cout << x << endl;

    Epetra_SerialDenseSolver solver;

    solver.SetMatrix(a);
    solver.SetVectors(x, b);
    solver.SolveWithTranspose(true);
    if (solver.ShouldEquilibrate())
    {
        solver.FactorWithEquilibration(true);
    }
    // solver.SolveToRefinedSolution(true);
    solver.Factor();
    solver.Solve();
    // solver.ApplyRefinement();
    solver.UnequilibrateLHS();
    
    if (!solver.Solved())// || !solver.SolutionRefined())
    {
        AssertMsg(false, "epetra failed to solve");
    }

    x_data.assign(solver.X(), solver.X() + size_);
}

// Direct matrix solution: faster than epetra for n<170
void Trilinos_Dense_Solve::
amesos_solve(vector<double> &a_data,
             vector<double> &b_data,
             vector<double> &x_data) const
{
    Check(a_data.size() == size_ * size_);
    Check(b_data.size() == size_);
    Check(x_data.size() == size_);

    int const index_base = 0;
    int const num_elements = size_;
    string solver_type = "Klu";
    
    vector<int> const num_entries_per_row(size_, size_);
    vector<int> column_indices(size_);
    for (unsigned i=0; i<size_; ++i)
    {
        column_indices[i] = i;
    }
    
#ifdef EPETRA_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif
    Epetra_Map map(num_elements, index_base, comm);
    Epetra_Vector lhs(View, map, &x_data[0]);
    Epetra_Vector rhs(View, map, &b_data[0]);
    Epetra_CrsMatrix matrix(View, map, size_);
    
    for (unsigned i=0; i<size_; ++i)
    {
        matrix.InsertGlobalValues(i, num_entries_per_row[i], &a_data[size_*i], &column_indices[0]);
    }
    matrix.FillComplete();
    
    Epetra_LinearProblem problem(&matrix, &lhs, &rhs);
    Amesos factory;
    
    Amesos_BaseSolver* solver = factory.Create(solver_type, problem);
    if (solver == NULL)
    {
        AssertMsg(false, "specified solver is not available");
    }
    
    solver->SymbolicFactorization();
    solver->NumericFactorization();
    solver->Solve();
}

// Does not work consistently for dense matrices
void Trilinos_Dense_Solve::
aztec_solve(vector<double> &a_data,
            vector<double> &b_data,
            vector<double> &x_data) const
{
    Check(a_data.size() == size_ * size_);
    Check(b_data.size() == size_);
    Check(x_data.size() == size_);

    int const index_base = 0;
    int const num_elements = size_;
    int const max_iterations = 10000;
    double const tolerance = 1.0E-6;
    
    vector<int> const num_entries_per_row(size_, size_);
    vector<int> column_indices(size_);
    for (unsigned i=0; i<size_; ++i)
    {
        column_indices[i] = i;
    }
    
#ifdef EPETRA_MPI
    Epetra_MpiComm comm( MPI_COMM_WORLD );
#else
    Epetra_SerialComm comm;
#endif
    Epetra_Map map(num_elements, index_base, comm);
    Epetra_Vector lhs(View, map, &x_data[0]);
    Epetra_Vector rhs(View, map, &b_data[0]);
    Epetra_CrsMatrix matrix(View, map, size_);

    for (unsigned i=0; i<size_; ++i)
    {
        matrix.InsertGlobalValues(i, num_entries_per_row[i], &a_data[size_*i], &column_indices[0]);
    }
    matrix.FillComplete();
    
    Epetra_LinearProblem problem(&matrix, &lhs, &rhs);
    AztecOO solver(problem);
    solver.SetAztecOption(AZ_precond, AZ_none);
    // solver.SetAztecOption(AZ_precond, AZ_Jacobi);
    // solver.SetAztecOption(AZ_poly_ord, 3);
    solver.SetAztecOption(AZ_solver, AZ_gmres);
    solver.SetAztecOption(AZ_kspace, 1000);
    solver.SetAztecOption(AZ_output, AZ_none);
    solver.Iterate(max_iterations, tolerance);
}
