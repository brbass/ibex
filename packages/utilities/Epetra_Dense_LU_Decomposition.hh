#ifndef Epetra_LU_Decomposition_hh
#define Epetra_LU_Decomposition_hh

#include "LU_Decomposition.hh"

class Epetra_SerialDenseMatrix;
class Epetra_SerialDenseSolver;

class Epetra_LU_Decomposition : public LU_Decomposition<double>
{
public:

    Epetra_LU_Decomposition(int size);
    
    // Check whether data has been initialized
    virtual bool initialized() const
    {
        return initialized_;
    }
    
    // Set matrix and perform decomposition
    virtual void initialize(std::vector<double> &a) const;
    
    // Rank of matrix
    virtual int size() const override
    {
        return size_;
    }

    // Apply to one vector
    virtual void solve(std::vector<double> &b,
                       std::vector<double> &x) const override;

    // Apply to multiple vectors (possibly a matrix)
    virtual void multi_solve(int number_of_vectors,
                             std::vector<double> &b,
                             std::vector<double> &x) const override;

    // Inverse of matrix
    virtual void inverse(std::vector<double> &ainv) const override;

    // Determinant of matrix
    virtual double determinant() const override;

private:

    bool initialized_;
    int size_;
    shared_ptr<Epetra_SerialDenseMatrix> a_;
    shared_ptr<Epetra_SerialDenseSolver> solver_;
    
};

#endif
