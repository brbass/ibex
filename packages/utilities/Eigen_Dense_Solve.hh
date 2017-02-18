#ifndef Eigen_Dense_Solve_hh
#define Eigen_Dense_Solve_hh

#include <vector>

#include <Eigen/Dense>

template<class Scalar>
class Eigen_Dense_Solve
{
public:
    
    typedef Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > EMatrixR;
    typedef Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > EVector;
    
    Eigen_Dense_Solve(int size):
        size_(size)
    {
    }
    
    virtual int size() const override
    {
        return size_;
    }
    
    void solve(std::vector<Scalar> &a_data,
               std::vector<Scalar> &b_data,
               std::vector<Scalar> &x_data)
    {
        Check(a_data.size() == size_ * size_);
        Check(b_data.size() == size_);
        Check(x_data.size() == size_);
        
        EMatrixR A(&a_data[0], size_, size_);
        EVector b(&b_data[0], size_);
        EVector x(&x_data[0], size_);
        
        x = A.fullPivLu().solve(b);
    }
};

#endif
