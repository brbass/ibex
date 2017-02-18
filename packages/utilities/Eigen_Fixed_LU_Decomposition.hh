#ifndef Eigen_Fixed_LU_Decomposition_hh
#define Eigen_Fixed_LU_Decomposition_hh

#include "LU_Decomposition.hh"

template<int const size_, class Scalar>
class Eigen_Fixed_LU_Decomposition
{
public:
    
    typedef Eigen::Map<Eigen::Matrix<Scalar, size_, size_, Eigen::RowMajor> const> EMatrixRC;
    typedef Eigen::Map<Eigen::Matrix<Scalar, size_, 1> const> EVectorC;
    typedef Eigen::Map<Eigen::Matrix<Scalar, size_, 1> > EVector;
    
    
};

#endif
