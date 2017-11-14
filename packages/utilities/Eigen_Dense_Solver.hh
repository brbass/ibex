#ifndef Eigen_Dense_Solver_hh
#define Eigen_Dense_Solver_hh

#include <vector>

#include <Eigen/Dense>

#include "Dense_Solver.hh"

template<class Scalar>
class Eigen_Dense_Solver : public Dense_Solver<Scalar>
{
public:
    
    // Matrices and vectors
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> EMatrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EMatrixR;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> EVector;
    
    // Mapped types
    typedef Eigen::Map<EVector> EMVector;
    typedef Eigen::Map<EMatrix> EMMatrix;
    typedef Eigen::Map<EMatrixR> EMMatrixR;

    // LU decomposition
    typedef Eigen::FullPivLU<EMatrixR> ELU;

    // Constructor
    Eigen_Dense_Solver(int size):
        initialized_(false),
        size_(size)
    {
    }

    // Rank of matrix
    virtual int size() const override
    {
        return size_;
    }
    
    // Check whether data has been initialized
    virtual bool initialized() const override
    {
        return initialized_;
    }
    
    // Set matrix and perform decomposition
    virtual void initialize(std::vector<Scalar> &a_data) override
    {
        Check(a_data.size() == size_ * size_);
        
        a_ = EMMatrixR(&a_data[0], size_, size_);
        lu_ = a_.fullPivLu();
        initialized_ = true;
    }
    
    // Solve problem Ax=b using temporary data (no initialization needed)
    virtual void solve(std::vector<Scalar> &a_data,
                       std::vector<Scalar> &b_data,
                       std::vector<Scalar> &x_data) override
    {
        Check(a_data.size() == size_ * size_);
        Check(b_data.size() == size_);
        Check(x_data.size() == size_);
        
        EMMatrixR a(&a_data[0], size_, size_);
        EMVector b(&b_data[0], size_);
        EMVector x(&x_data[0], size_);
        
        x = a.fullPivLu().solve(b);
    }
    
    // Apply to one vector
    virtual void solve(std::vector<double> &b_data,
                       std::vector<double> &x_data) override
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
                             std::vector<double> &b_data,
                             std::vector<double> &x_data) override
    {
        Assert(initialized_);
        Check(b_data.size() == size_ * number_of_vectors);
        Check(x_data.size() == size_ * number_of_vectors);

        EMMatrix b(&b_data[0], size_, number_of_vectors);
        EMMatrix x(&x_data[0], size_, number_of_vectors);

        x = lu_.solve(b);
    }

    // Get inverse
    virtual void inverse(std::vector<double> &ainv_data) override
    {
        Assert(initialized_);
        Check(ainv_data.size() == size_ * size_);
        
        EMMatrixR ainv(&ainv_data[0], size_, size_);

        ainv = lu_.inverse();
    }

    virtual void inverse(std::vector<double> &a_data,
                         std::vector<double> &ainv_data) override
    {
        Check(a_data.size() == size_ * size_);
        Check(ainv_data.size() == size_ * size_);

        EMMatrixR a(&a_data[0], size_, size_);
        EMMatrixR ainv(&ainv_data[0], size_, size_);

        ainv = a.fullPivLu().inverse();
    }
    
    // Get determinant
    virtual Scalar determinant() override
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
