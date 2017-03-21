#ifndef Vector_Operator_Difference_hh
#define Vector_Operator_Difference_hh

#include "Vector_Operator.hh"

#include <memory>

/* 
   Gives the difference of two vector operators,
   
   op1(x) - op2(x)
*/
class Vector_Operator_Difference : public Vector_Operator
{
public:

    Vector_Operator_Difference(std::shared_ptr<Vector_Operator> op1,
                               std::shared_ptr<Vector_Operator> op2);

    virtual int row_size() const override
    {
        return op1_->row_size();
    }
    virtual int column_size() const override
    {
        return op1_->column_size();
    }
    
    virtual void check_class_invariants() const override;
    
private:

    virtual void apply(std::vector<double> &x) const override;

    std::shared_ptr<Vector_Operator> op1_;
    std::shared_ptr<Vector_Operator> op2_;
};

#endif
