#ifndef Dense_Solve_hh
#define Dense_Solve_hh

#include <memory>
#include <vector>

/*
  Solve the problem Ax=b, discarding any intermediate results
*/
template<class Scalar>
class Dense_Solve
{
public:

    // Creator
    Dense_Solve()
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
