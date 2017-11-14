#include "Epetra_Dense_Solver.hh"

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseSolver.h>

#include "Check.hh"

using namespace std;

Epetra_Dense_Solver::
Epetra_Dense_Solver(int size):
    initialized_(false),
    size_(size)
{
}

void Epetra_Dense_Solver::
initialize(vector<double> &a_data)
{
    Check(a_data.size() == size_ * size_);
    
    a_ = make_shared<Epetra_SerialDenseMatrix>(Copy,
                                               &a_data[0],
                                               size_,
                                               size_,
                                               size_);
    solver_ = make_shared<Epetra_SerialDenseSolver>();
    solver_->SetMatrix(*a_);
    solver_->SolveWithTranspose(true);
    solver_->Factor();
    initialized_ = true;
}

void Epetra_Dense_Solver::
solve(vector<double> &a_data,
      vector<double> &b_data,
      vector<double> &x_data)
{
    Check(a_data.size() == size_ * size_);
    Check(b_data.size() == size_);
    Check(x_data.size() == size_);

    Epetra_SerialDenseMatrix a(View, &a_data[0], size_, size_, size_);
    Epetra_SerialDenseVector b(View, &b_data[0], size_);
    Epetra_SerialDenseVector x(size_);

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

void Epetra_Dense_Solver::
solve(vector<double> &b_data,
      vector<double> &x_data)
{
    Assert(initialized_);
    Check(b_data.size() == size_);
    Check(x_data.size() == size_);
    
    Epetra_SerialDenseVector b(View, &b_data[0], size_);
    Epetra_SerialDenseVector x(size_);
    solver_->SetVectors(x, b);
    solver_->Solve();

    if (!solver_->Solved())
    {
        AssertMsg(false, "epetra failed to solve");
    }
    
    x_data.assign(solver_->X(), solver_->X() + size_);
}

void Epetra_Dense_Solver::
multi_solve(int number_of_vectors,
            std::vector<double> &b_data,
            std::vector<double> &x_data)
{
    Assert(initialized_);
    Check(b_data.size() == size_ * number_of_vectors);
    Check(x_data.size() == size_ * number_of_vectors);

    Epetra_SerialDenseMatrix b(View,
                               &b_data[0],
                               size_,
                               size_,
                               number_of_vectors);
    Epetra_SerialDenseMatrix x(size_,
                               number_of_vectors);
    solver_->SetVectors(x, b);
    solver_->Solve();
    
    if (!solver_->Solved())
    {
        AssertMsg(false, "epetra failed to solve");
    }
    
    x_data.assign(solver_->X(), solver_->X() + size_ * number_of_vectors);
}

void Epetra_Dense_Solver::
inverse(std::vector<double> &ainv)
{
    Assert(initialized_);
    Assert(ainv.size() == size_ * size_);
    AssertMsg(false, "not implemented");
}

void Epetra_Dense_Solver::
inverse(std::vector<double> &a,
        std::vector<double> &ainv)
{
    Assert(initialized_);
    Assert(a.size() == size_ * size_);
    Assert(ainv.size() == size_ * size_);
    AssertMsg(false, "not implemented");
}

double Epetra_Dense_Solver::
determinant()
{
    Assert(initialized_);
    AssertMsg(false, "not implemented");
    return -1;
}
