#ifndef Epetra_Dense_Solver_hh
#define Epetra_Dense_Solver_hh

#include <memory>

#include "Dense_Solver.hh"

class Epetra_SerialDenseMatrix;
class Epetra_SerialDenseSolver;

/*
  Epetra dense solver
  Aztec (iterative) and Amesos (sparse) performed poorly for large, dense problems and were removed
*/
class Epetra_Dense_Solver : public Dense_Solver<double>
{
public:

    // Constructor
    Epetra_Dense_Solver(int size);
    

    // Check whether data has been initialized
    virtual bool initialized() const override
    {
        return initialized_;
    }
    
    // Set matrix and perform decomposition
    virtual void initialize(std::vector<double> &a) override;
    
    // Rank of matrix
    virtual int size() const override
    {
        return size_;
    }
 
    // Solve problem Ax=b using temporary data (no initialization needed)
    virtual void solve(std::vector<double> &a_data,
                       std::vector<double> &b_data,
                       std::vector<double> &x_data) override;
    
   // Apply to one vector
    virtual void solve(std::vector<double> &b,
                       std::vector<double> &x) override;

    // Apply to multiple vectors (possibly a matrix)
    virtual void multi_solve(int number_of_vectors,
                             std::vector<double> &b,
                             std::vector<double> &x) override;

    // Inverse of matrix
    virtual void inverse(std::vector<double> &ainv) override;

    // Determinant of matrix
    virtual double determinant() override;

private:

    bool initialized_;
    int size_;
    std::shared_ptr<Epetra_SerialDenseMatrix> a_;
    std::shared_ptr<Epetra_SerialDenseSolver> solver_;
};

#endif
