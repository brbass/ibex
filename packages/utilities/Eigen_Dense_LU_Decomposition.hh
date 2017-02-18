#ifndef Eigen_Dense_LU_Decomposition_hh
#define Eigen_Dense_LU_Decomposition_hh

#include <vector>

#include <Eigen/Dense>

#include "LU_Decomposition.hh"

template<class Scalar>
class Eigen_Dense_LU_Decomposition : public LU_Decomposition<Scalar>
{
public:

    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EMatrixR;
    typedef Eigen::Map<EMatrixR> EMMatrixR;
    typedef Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > EMVector;
    typedef Eigen::FullPivLU<EMatrixR> ELU;
    
    Eigen_Dense_LU_Decomposition(int size):
        initialized_(false),
        size_(size)
    {
    }

    // Check whether data has been initialized
    virtual bool initialized() const override
    {
        return initialized_;
    }
    
    // Set matrix and perform decomposition
    virtual void initialize(std::vector<Scalar> &a) const override
    {
        Check(a_data.size() == size_ * size_);
        
        a_ = EMMatrixR(&a_data[0], size_, size_);
        lu_ = a_.fullPivLu();
        initialized_ = true;
    }
    
    // Rank of matrix
    virtual int size() const override
    {
        return size_;
    }

    // Apply to one vector
    virtual void solve(std::vector<double> &b_data,
                       std::vector<double> &x_data) const override
    {
        Assert(initialized_);
        Check(b_data.size() == size_);
        Check(x_data.size() == size_);
        
        EMVector b(&b_data[0], size_);
        EMVector x(&x_data[0], size_);
        
        x = lu_.solve(b);
    }
    
    // Apply to multiple vectors (possibly a matrix)
    virtual void multi_solve(int number_of_vectors,
                             std::vector<double> &b,
                             std::vector<double> &x) const override
    {
        Assert(initialized_);
        AssertMsg(false, "not implemented");
    }

    // Get inverse
    virtual void inverse(std::vector<double> &ainv_data) const override
    {
        Assert(initialized_);
        Check(ainv_data.size() == size_ * size_);
        
        EMMatrixR ainv(&ainv_data[0], size_, size_);

        ainv = lu_.inverse();
    }

    // Get determinant
    virtual Scalar determinant() const override
    {
        Assert(initialized_);
        
        return lu_.determinant();
    }
    
private:

    bool initialized_;
    int size_;
    EMatrixR a_;
    ELU lu_;
};

#endif
