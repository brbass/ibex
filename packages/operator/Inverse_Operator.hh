#ifndef Inverse_Operator_hh
#define Inverse_Operator_hh

#include <memory>

#include "Square_Vector_Operator.hh"

/*
  Inverts a Vector_Operator
*/
class Inverse_Operator : public Square_Vector_Operator
{
public:
    
    // Constructor
    Inverse_Operator(std::shared_ptr<Vector_Operator> vector_operator);
    
    virtual void check_class_invariants() const override = 0;
    
    virtual int size() const override
    {
        return size_;
    }
    virtual std::string description() const override = 0;
    
protected:

    int size_;
    std::shared_ptr<Vector_Operator> vector_operator_;
    
    
private:
    
    virtual void apply(std::vector<double> &x) const override = 0;
};

#endif
