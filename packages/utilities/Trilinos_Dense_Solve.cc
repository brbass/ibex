#include "Trilinos_Dense_Solve.hh"

#include <string>
#include <vector>

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseSolver.h>

#include <Amesos.h>
#include <AztecOO.h>
// #include <AztecOO_Version.h>
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
Trilinos_Dense_Solve()
{
}

// Iterative matrix solution: faster than GSL for n>300
void Trilinos_Dense_Solve::
epetra_solve(vector<double> &a_data,
             vector<double> &b_data,
             vector<double> &x_data,
             unsigned number_of_elements)
{
    Check(a_data.size() == number_of_elements*number_of_elements);
    Check(b_data.size() == number_of_elements);
    Check(x_data.size() == number_of_elements);

    Epetra_SerialDenseMatrix a(View, &a_data[0], number_of_elements, number_of_elements, number_of_elements);
    Epetra_SerialDenseVector x(number_of_elements);
    Epetra_SerialDenseVector b(View, &b_data[0], number_of_elements);

    
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

    x_data.assign(solver.X(), solver.X() + number_of_elements);
}

// Direct matrix solution: faster than epetra for n<170
void Trilinos_Dense_Solve::
amesos_dense_solve(vector<double> &a_data,
                   vector<double> &b_data,
                   vector<double> &x_data,
                   unsigned number_of_elements)
{
    Check(a_data.size() == number_of_elements*number_of_elements);
    Check(b_data.size() == number_of_elements);
    Check(x_data.size() == number_of_elements);

    int const index_base = 0;
    int const num_elements = number_of_elements;
    string solver_type = "Klu";
    
    vector<int> const num_entries_per_row(number_of_elements, number_of_elements);
    vector<int> column_indices(number_of_elements);
    for (unsigned i=0; i<number_of_elements; ++i)
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
    Epetra_CrsMatrix matrix(View, map, number_of_elements);

    for (unsigned i=0; i<number_of_elements; ++i)
    {
        matrix.InsertGlobalValues(i, num_entries_per_row[i], &a_data[number_of_elements*i], &column_indices[0]);
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
aztec_dense_solve(vector<double> &a_data,
                  vector<double> &b_data,
                  vector<double> &x_data,
                  unsigned number_of_elements)
{
    Check(a_data.size() == number_of_elements*number_of_elements);
    Check(b_data.size() == number_of_elements);
    Check(x_data.size() == number_of_elements);

    int const index_base = 0;
    int const num_elements = number_of_elements;
    int const max_iterations = 10000;
    double const tolerance = 1.0E-6;
    
    vector<int> const num_entries_per_row(number_of_elements, number_of_elements);
    vector<int> column_indices(number_of_elements);
    for (unsigned i=0; i<number_of_elements; ++i)
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
    Epetra_CrsMatrix matrix(View, map, number_of_elements);

    for (unsigned i=0; i<number_of_elements; ++i)
    {
        matrix.InsertGlobalValues(i, num_entries_per_row[i], &a_data[number_of_elements*i], &column_indices[0]);
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
