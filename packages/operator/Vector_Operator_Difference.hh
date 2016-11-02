#ifndef Vector_Operator_Difference_hh
#define Vector_Operator_Difference_hh

#include "Vector_Operator.hh"

#include <memory>

using std::shared_ptr;

/* 
   Gives the difference of two vector operators,
   
   op1(x) - op2(x)
*/
class Vector_Operator_Difference : public Vector_Operator
{
public:

    Vector_Operator_Difference(shared_ptr<Vector_Operator> op1,
                               shared_ptr<Vector_Operator> op2);

    virtual void check_class_invariants() const override;
    
private:

    virtual void apply(vector<double> &x) const override;

    shared_ptr<Vector_Operator> op1_;
    shared_ptr<Vector_Operator> op2_;
};

#endif
