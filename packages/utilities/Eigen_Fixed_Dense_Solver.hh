#ifndef Eigen_Fixed_Dense_Solver_hh
#define Eigen_Fixed_Dense_Solver_hh

#include <vector>

#include <Eigen/Dense>

#include "Dense_Solver.hh"

template<int const size_, class Scalar>
class Eigen_Fixed_Dense_Solver : public Dense_Solver<Scalar>
{
public:
    
    typedef Eigen::Map<Eigen::Matrix<Scalar, size_, size_, Eigen::RowMajor> > EMatrixR;
    typedef Eigen::Map<Eigen::Matrix<Scalar, size_, 1> > EVector;
    
    Eigen_Fixed_Dense_Solver()
    {
    }

    virtual int size() const override
    {
        return size_;
    }
    
    virtual void solve(std::vector<Scalar> &a_data,
                       std::vector<Scalar> &b_data,
                       std::vector<Scalar> &x_data) const override
    {
        Check(a_data.size() == size_ * size_);
        Check(b_data.size() == size_);
        Check(x_data.size() == size_);
        
        EMatrixR a(&a_data[0]);
        EVector b(&b_data[0]);
        EVector x(&x_data[0]);
        
        x = a.fullPivLu().solve(b);
    }
};

#endif
