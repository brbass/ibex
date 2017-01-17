#ifndef Vector_Operator_Product_hh
#define Vector_Operator_Product_hh

#include "Vector_Operator.hh"

#include <memory>

/* 
   Gives the product of two vector operators,
   
   op1(op2(x))
*/
class Vector_Operator_Product : public Vector_Operator
{
public:

    Vector_Operator_Product(std::shared_ptr<Vector_Operator> op1,
                            std::shared_ptr<Vector_Operator> op2);

    virtual void check_class_invariants() const override;
    
private:

    virtual void apply(vector<double> &x) const override;

    std::shared_ptr<Vector_Operator> op1_;
    std::shared_ptr<Vector_Operator> op2_;
};

#endif
