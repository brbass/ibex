#ifndef Aztec_Inverse_Operator_hh
#define Aztec_Inverse_Operator_hh

#include <memory>

#include "Inverse_Operator.hh"

class AztecOO;
class Epetra_Comm;
class Epetra_LinearProblem;
class Epetra_Map;
class Epetra_Operator_Interface;
class Epetra_Vector;

/*
  Inverts a Vector_Operator
*/
class Aztec_Inverse_Operator : public Inverse_Operator
{
public:

    struct Options
    {
        int max_iterations = 1000;
        int kspace = 20;
        int solver_print = 0;
        double tolerance = 1e-10;
    };
    
    // Constructor
    Aztec_Inverse_Operator(Options options,
                           std::shared_ptr<Vector_Operator> vector_operator);
    Aztec_Inverse_Operator(Options options,
                           std::shared_ptr<Vector_Operator> vector_operator,
                           std::shared_ptr<Epetra_Comm> comm,
                           std::shared_ptr<Epetra_Map> map);
    
    
    virtual void check_class_invariants() const override;
    
    virtual std::string description() const override;
    
private:
    
    virtual void apply(std::vector<double> &x) const override;
    virtual void initialize_trilinos(bool initialize_comm);
    
    // Solver data
    Options options_;
    std::shared_ptr<Epetra_Comm> comm_;
    std::shared_ptr<Epetra_Map> map_;
    std::shared_ptr<Epetra_Vector> lhs_;
    std::shared_ptr<Epetra_Vector> rhs_;
    std::shared_ptr<Epetra_Operator_Interface> oper_;
    std::shared_ptr<Epetra_LinearProblem> problem_;
    std::shared_ptr<AztecOO> solver_;
};

#endif
