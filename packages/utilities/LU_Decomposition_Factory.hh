#ifndef LU_Decomposition_Factory_hh
#define LU_Decomposition_Factory_hh

#include <memory>

#include "LU_Decomposition.hh"

class LU_Decomposition_Factory
{
public:

    typedef double Scalar;
    
    enum class Type
    {
        DEFAULT,
        EIGEN,
        EIGEN_FIXED,
        EPETRA
    };
    
    LU_Decomposition_Factory();
    
    std::shared_ptr<LU_Decomposition<Scalar> > get_solver(int size,
                                                          Type type = Type::DEFAULT) const;

private:
    
    std::shared_ptr<LU_Decomposition<Scalar> > get_default(int size) const;
    std::shared_ptr<LU_Decomposition<Scalar> > get_fixed_eigen_solver(int size) const;

};

#endif
