#ifndef Eigen_Dense_Solver_hh
#define Eigen_Dense_Solver_hh

#include <vector>

#include <Eigen/Dense>

#include "Dense_Solver.hh"

template<class Scalar>
class Eigen_Dense_Solver : public Dense_Solver<Scalar>
{
public:
    
    typedef Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > EMatrixR;
    typedef Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > EVector;
    
    Eigen_Dense_Solver(int size):
        size_(size)
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
        
        EMatrixR a(&a_data[0], size_, size_);
        EVector b(&b_data[0], size_);
        EVector x(&x_data[0], size_);
        
        x = a.fullPivLu().solve(b);
    }

private:

    int size_;
};

#endif
