#ifndef Eigen_Fixed_Dense_Solve_hh
#define Eigen_Fixed_Dense_Solve_hh

#include <vector>

#include <Eigen/Dense>

template<int const size_, class Scalar>
class Eigen_Fixed_Dense_Solve
{
public:
    
    typedef Eigen::Map<Eigen::Matrix<Scalar, size_, size_, Eigen::RowMajor> const> EMatrixRC;
    typedef Eigen::Map<Eigen::Matrix<Scalar, size_, 1> const> EVectorC;
    typedef Eigen::Map<Eigen::Matrix<Scalar, size_, 1> > EVector;
    
    Eigen_Dense_Solve()
    {
    }

    virtual int size() const override
    {
        return size_;
    }
    
    virtual void solve(std::vector<Scalar> const &a_data,
                       std::vector<Scalar> const &b_data,
                       std::vector<Scalar> &x_data) const override
    {
        Check(a_data.size() == size_ * size_);
        Check(b_data.size() == size_);
        Check(x_data.size() == size_);
        
        EMatrixRC A(&a_data[0]);
        EVectorC b(&b_data[0]);
        EVector x(&x_data[0]);
        
        x = A.fullPivLu().solve(b);
    }
};

