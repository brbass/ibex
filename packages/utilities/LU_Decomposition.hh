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

    // Rank of matrix
    virtual int size() const = 0;

    // Apply to one vector
    virtual void apply(std::vector<double> &b,
                       std::vector<double> &x) const = 0;

    // Apply to multiple vectors (possibly a matrix)
    virtual void multi_apply(int number_of_vectors,
                             std::vector<double> &b,
                             std::vector<double> &x) const = 0;
};

#endif
