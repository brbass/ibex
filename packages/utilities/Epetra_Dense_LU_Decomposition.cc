#include "Epetra_Dense_LU_Decomposition.hh"

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseSolver.h>

#include "Check.hh"

Epetra_LU_Decomposition::
Epetra_LU_Decomposition(int size):
    initialized_(false),
    size_(size)
{
}


void Epetra_LU_Decomposition::
initialize(vector<double> &a_data)
{
    Check(a_data.size() == size_ * size_);
    
    a_ = make_shared<Epetra_SerialDenseMatrix>(Copy,
                                               &a_data[0],
                                               size_,
                                               size_,
                                               size_);
    solver_ = make_shared<Epetra_SerialDenseSolver>();
    solver_->SetMatrix(a_);
    solver_->SolveWithTranspose(true);
    solver_->Factor();
    initialized_ = true;
}

void Epetra_LU_Decomposition::
solve(vector<double> &b_data,
      vector<double> &x_data) const
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
    
    x_data.assign(solver.X(), solver.X() + size_);
}

void Epetra_LU_Decomposition::
multi_solve(int number_of_vectors,
            std::vector<double> &b,
            std::vector<double> &x)
{
    Assert(initialized_);
    Check(b_data.size() == size_ * number_of_vectors);
    Check(x_data.size() == size_ * number_of_vectors);
    AssertMsg(false, "not implemented");
}

void Epetra_LU_Decomposition::
inverse(std::vector<double> &ainv) const
{
    Assert(initialized_);
    Assert(ainv.size() == size_ * size_);
    AssertMsg(false, "not implemented");
}

double Epetra_LU_Decomposition::
determinant() const
{
    Assert(initialized_);
    AssertMsg(false, "not implemented");
}
