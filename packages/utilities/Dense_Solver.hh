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

    // Creator
    Dense_Solver()
    {
    }

    // Solve problem Ax=b
    virtual void solve(std::vector<Scalar> &a_data, 
                       std::vector<Scalar> &b_data,
                       std::vector<Scalar> &x_data) const = 0;
    
    // Size of square matrix
    virtual int size() const = 0;
};

#endif
