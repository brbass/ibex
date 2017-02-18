#ifndef Direct_Dense_Solver_hh
#define Direct_Dense_Solver_hh

#include <string>
#include <vector>

#include "Dense_Solver.hh"
#include "Linear_Algebra.hh"

template<class Scalar>
class Direct_Dense_Solver : public Dense_Solver<Scalar>
{
public:
    
    Direct_Dense_Solver(int size):
        size_(size)
    {
    }
    
    virtual void solve(std::vector<Scalar> &a_data,
                       std::vector<Scalar> &b_data,
                       std::vector<Scalar> &x_data) const override
    {
        Check(a_data.size() == size_ * size_);
        Check(b_data.size() == size_);
        Check(x_data.size() == size_);
        
        switch(size_)
        {
        case 1:
            return Linear_Algebra::linear_solve_1(a_data, b_data, x_data);
        case 2:
            return Linear_Algebra::linear_solve_2(a_data, b_data, x_data);
        case 3:
            return Linear_Algebra::linear_solve_3(a_data, b_data, x_data);
        case 4:
            return Linear_Algebra::linear_solve_4(a_data, b_data, x_data);
        case 5:
            return Linear_Algebra::linear_solve_5(a_data, b_data, x_data);
        default:
            AssertMsg(false, "linear solve of size (" + std::to_string(size_) + ") not implemented");
            return;
        }
    }
    
    virtual int size() const override
    {
        return size_;
    }

private:

    int size_;
};

#endif
