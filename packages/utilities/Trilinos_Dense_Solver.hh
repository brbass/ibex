#ifndef Trilinos_Dense_Solver_hh
#define Trilinos_Dense_Solver_hh

#include "Dense_Solver.hh"

/*
  Epetra dense solver
  Aztec (iterative) and Amesos (sparse) performed poorly for large, dense problems and were removed
*/
class Trilinos_Dense_Solver : public Dense_Solver<double>
{
public:

    Trilinos_Dense_Solver(int size);
    
    // Solve using Epetra_SerialDenseSolver (LU decomposition)
    virtual void solve(std::vector<double> &a_data,
                       std::vector<double> &b_data,
                       std::vector<double> &x_data) const override;

    virtual int size() const override
    {
        return size_;
    }
    
private:

    int size_;
};

#endif
