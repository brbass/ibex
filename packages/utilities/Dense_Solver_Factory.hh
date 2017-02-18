#ifndef Dense_Solver_Factory_hh
#define Dense_Solver_Factory_hh

#include <memory>

#include "Dense_Solver.hh"

class Dense_Solver_Factory
{
public:

    typedef double Scalar;
    
    enum class Type
    {
        DEFAULT,
        DIRECT,
        EIGEN,
        EIGEN_FIXED,
        EPETRA
    };
    
    Dense_Solver_Factory();
    
    std::shared_ptr<Dense_Solver<Scalar> > get_solver(int size,
                                                      Type type = Type::DEFAULT) const;

private:
    
    std::shared_ptr<Dense_Solver<Scalar> > get_default(int size) const;
    std::shared_ptr<Dense_Solver<Scalar> > get_fixed_eigen_solver(int size) const;

};

#endif
