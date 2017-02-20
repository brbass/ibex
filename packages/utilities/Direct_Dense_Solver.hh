#ifndef Direct_Dense_Solver_hh
#define Direct_Dense_Solver_hh

#include <string>
#include <vector>

#include "Dense_Solver.hh"
#include "Linear_Algebra.hh"

/*
  Class to allow partial template class specialization
  Contains shared functions
*/
template<int const size_, class Scalar>
class Direct_Dense_Solver_Implementation : public Dense_Solver<Scalar>
{
public:

    // Constructor
    Direct_Dense_Solver_Implementation()
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
        Check(a.size() == size_ * size_);
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

        return solve_(a_data,
                      b_data,
                      x_data);
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
        return inverse_(ainv_data);
    }

    // Determinant of matrix
    virtual Scalar determinant() override
    {
        return determinant_();
    }
    
protected:

    // Data
    bool initialized_;
    std::vector<Scalar> a_;

private:
    
    // Implement these functions for each desired size
    virtual void solve_(std::vector<Scalar> &a_data,
                        std::vector<Scalar> &b_data,
                        std::vector<Scalar> &x_data) = 0;
    virtual void inverse_(std::vector<Scalar> &ainv_data) = 0;
    virtual Scalar determinant_() = 0;
    
};

/*
  Template for direct dense solver
  Must be specialized for each desired size
*/
template<int const size_, class Scalar>
class Direct_Dense_Solver : public Direct_Dense_Solver_Implementation<size_, Scalar>
{
private:
    virtual void solve_(std::vector<Scalar> &a_data,
                        std::vector<Scalar> &b_data,
                        std::vector<Scalar> &x_data) override
    {
        AssertMsg(false, "not implemented for size " + std::to_string(size_));
    }
    virtual void inverse_(std::vector<Scalar> &ainv_data) override
    {
        AssertMsg(false, "not implemented for size " + std::to_string(size_));
    }
    virtual Scalar determinant_() override
    {
        AssertMsg(false, "not implemented for size " + std::to_string(size_));
    }
};

template<class Scalar>
class Direct_Dense_Solver<1, Scalar> : public Direct_Dense_Solver_Implementation<1, Scalar>
{
private:
    virtual void solve_(std::vector<Scalar> &a,
                        std::vector<Scalar> &b,
                        std::vector<Scalar> &x)
    {
        return Linear_Algebra::linear_solve_1(a,
                                              b,
                                              x);
    }
    virtual void inverse_(std::vector<Scalar> &ainv_data)
    {
        return Linear_Algebra::direct_inverse_1(this->a_,
                                                ainv_data);
    }

    virtual Scalar determinant_()
    {
        return Linear_Algebra::determinant_1(this->a_);
    }
};

template<class Scalar>
class Direct_Dense_Solver<2, Scalar> : public Direct_Dense_Solver_Implementation<2, Scalar>
{
private:
    virtual void solve_(std::vector<Scalar> &a,
                        std::vector<Scalar> &b,
                        std::vector<Scalar> &x)
    {
        return Linear_Algebra::linear_solve_2(a,
                                              b,
                                              x);
    }
    virtual void inverse_(std::vector<Scalar> &ainv_data)
    {
        return Linear_Algebra::direct_inverse_2(this->a_,
                                                ainv_data);
    }

    virtual Scalar determinant_()
    {
        return Linear_Algebra::determinant_2(this->a_);
    }
};

template<class Scalar>
class Direct_Dense_Solver<3, Scalar> : public Direct_Dense_Solver_Implementation<3, Scalar>
{
private:
    virtual void solve_(std::vector<Scalar> &a,
                        std::vector<Scalar> &b,
                        std::vector<Scalar> &x)
    {
        return Linear_Algebra::linear_solve_3(a,
                                              b,
                                              x);
    }
    virtual void inverse_(std::vector<Scalar> &ainv_data)
    {
        return Linear_Algebra::direct_inverse_3(this->a_,
                                                ainv_data);
    }

    virtual Scalar determinant_()
    {
        return Linear_Algebra::determinant_3(this->a_);
    }
};

template<class Scalar>
class Direct_Dense_Solver<4, Scalar> : public Direct_Dense_Solver_Implementation<4, Scalar>
{
private:
    virtual void solve_(std::vector<Scalar> &a,
                        std::vector<Scalar> &b,
                        std::vector<Scalar> &x)
    {
        return Linear_Algebra::linear_solve_4(a,
                                              b,
                                              x);
    }
    virtual void inverse_(std::vector<Scalar> &ainv_data)
    {
        return Linear_Algebra::direct_inverse_4(this->a_,
                                                ainv_data);
    }

    virtual Scalar determinant_()
    {
        return Linear_Algebra::determinant_4(this->a_);
    }
};

template<class Scalar>
class Direct_Dense_Solver<5, Scalar> : public Direct_Dense_Solver_Implementation<5, Scalar>
{
private:
    virtual void solve_(std::vector<Scalar> &a,
                        std::vector<Scalar> &b,
                        std::vector<Scalar> &x)
    {
        return Linear_Algebra::linear_solve_5(a,
                                              b,
                                              x);
    }
    virtual void inverse_(std::vector<Scalar> &ainv_data)
    {
        return Linear_Algebra::direct_inverse_5(this->a_,
                                                ainv_data);
    }

    virtual Scalar determinant_()
    {
        return Linear_Algebra::determinant_5(this->a_);
    }
};

#endif
