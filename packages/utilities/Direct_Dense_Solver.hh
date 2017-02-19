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
    
    // Check whether data has been initialized
    virtual bool initialized() const
    {
        return initialized_;
    }
    
    // Set matrix and perform decomposition
    virtual void initialize(std::vector<Scalar> &a)
    {
        a_ = a;
        initialized_ = true;
    }

    // Rank of matrix
    virtual int size() const override
    {
        return size_;
    }

    // Solve problem Ax=b using temporary data (no initialization needed)
    virtual void solve(std::vector<Scalar> &a_data,
                       std::vector<Scalar> &b_data,
                       std::vector<Scalar> &x_data) override
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

    // Apply to one vector
    virtual void solve(std::vector<Scalar> &b_data,
                       std::vector<Scalar> &x_data) override
    {
        return solve(a_,
                     b_data,
                     x_data);
    }

    // Apply to multiple vectors (possibly a matrix)
    virtual void multi_solve(int number_of_vectors,
                             std::vector<Scalar> &b_data,
                             std::vector<Scalar> &x_data) override
    {
        AssertMsg(false, "not implemented");
    }
    
    // Inverse of matrix
    virtual void inverse(std::vector<Scalar> &ainv_data) override
    {
        Check(ainv_data.size() == size_ * size_);
        
        switch(size_)
        {
        case 1:
            return Linear_Algebra::direct_inverse_1(ainv_data, a_);
        case 2:
            return Linear_Algebra::direct_inverse_2(ainv_data, a_);
        case 3:
            return Linear_Algebra::direct_inverse_3(ainv_data, a_);
        case 4:
            return Linear_Algebra::direct_inverse_4(ainv_data, a_);
        case 5:
            return Linear_Algebra::direct_inverse_5(ainv_data, a_);
        default:
            AssertMsg(false, "linear solve of size (" + std::to_string(size_) + ") not implemented");
            return;
        }
    }

    // Determinant of matrix
    virtual Scalar determinant() override
    {
        AssertMsg(false, "not implemented");
    }
    
private:

    bool initialized_;
    int size_;
    std::vector<Scalar> a_;
};

#endif
