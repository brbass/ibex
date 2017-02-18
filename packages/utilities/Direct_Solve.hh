#ifndef Direct_Solve_hh
#define Direct_Solve_hh

#include <vector>

#include "Linear_Algebra.hh"

template<class Scalar>
class Direct_Solve : public Dense_Solve<Scalar>
{
public:
    
    Direct_Solve(int size):
        size_(size)
    {
    }
    
    virtual void solve(std::vector<Scalar> &a_data,
                       std::vector<Scalar> &b_data,
                       std::vector<Scalar> &x_data) const override;
    {
        Check(a_data.size() == size_ * size_);
        Check(b_data.size() == size_);
        Check(x_data.size() == size_);
        
        switch(size_)
        {
        case 1:
            return Linear_Algebra::linear_solve_1(a, b, x);
        case 2:
            return Linear_Algebra::linear_solve_2(a, b, x);
        case 3:
            return Linear_Algebra::linear_solve_3(a, b, x);
        case 4:
            return Linear_Algebra::linear_solve_4(a, b, x);
        case 5:
            return Linear_Algebra::linear_solve_5(a, b, x);
        default:
            AssertMsg(false, "linear solve of that size not implemented");
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
