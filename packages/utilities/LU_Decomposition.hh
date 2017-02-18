#ifndef LU_Decomposition_hh
#define LU_Decomposition_hh

#include <vector>

/*
  Calculate LU decomposition and apply to one or more vectors
*/
template<class Scalar>
class LU_Decomposition
{
public:

    LU_Decomposition()
    {
    }

    // Check whether data has been initialized
    virtual bool initialized() const = 0;
    
    // Set matrix and perform decomposition
    virtual void initialize(std::vector<Scalar> &a) const = 0;
    
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
};

#endif
