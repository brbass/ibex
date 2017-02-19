#ifndef Dense_Solver_hh
#define Dense_Solver_hh

#include <memory>
#include <vector>

/*
  Solve the problem Ax=b, discarding any intermediate results
*/
template<class Scalar>
class Dense_Solver
{
public:

    // Constructor
    Dense_Solver()
    {
    }
    
    // Check whether data has been initialized
    virtual bool initialized() const = 0;
    
    // Set matrix and perform decomposition
    virtual void initialize(std::vector<Scalar> &a) = 0;
    
    // Rank of matrix
    virtual int size() const = 0;

    // Apply to one vector
    virtual void solve(std::vector<Scalar> &b,
                       std::vector<Scalar> &x) const = 0;

    // Apply to multiple vectors (possibly a matrix)
    virtual void multi_solve(int number_of_vectors,
                             std::vector<Scalar> &b,
                             std::vector<Scalar> &x) const = 0;

    // Inverse of matrix
    virtual void inverse(std::vector<Scalar> &ainv) const = 0;

    // Determinant of matrix
    virtual Scalar determinant() const = 0;
    
    // Solve problem Ax=b using temporary data (no initialization needed)
    virtual void solve(std::vector<Scalar> &a_data, 
                       std::vector<Scalar> &b_data,
                       std::vector<Scalar> &x_data) const = 0;
    
    // Size of square matrix
    virtual int size() const = 0;
};

#endif
