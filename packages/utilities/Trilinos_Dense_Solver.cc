#include "Trilinos_Dense_Solver.hh"

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

Trilinos_Dense_Solver::
Trilinos_Dense_Solver(int size):
    size_(size)
{
}

void Trilinos_Dense_Solver::
solve(vector<double> &a_data,
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
    // if (solver.ShouldEquilibrate())
    // {
    //     solver.FactorWithEquilibration(true);
    // }
    // solver.SolveToRefinedSolution(true);
    solver.Factor();
    solver.Solve();
    // solver.ApplyRefinement();
    // solver.UnequilibrateLHS();
    
    if (!solver.Solved())// || !solver.SolutionRefined())
    {
        AssertMsg(false, "epetra failed to solve");
    }

    x_data.assign(solver.X(), solver.X() + size_);
}

