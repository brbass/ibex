#include "Vector_Operator_Functions.hh"

#include "Vector_Operator_Difference.hh"
#include "Vector_Operator_Product.hh"
#include "Vector_Operator_Sum.hh"

using std::make_shared;
using std::shared_ptr;

shared_ptr<Vector_Operator> operator+(shared_ptr<Vector_Operator> const op1,
                                      shared_ptr<Vector_Operator> const op2)
{
    return make_shared<Vector_Operator_Sum>(op1, op2);
}

// Subtract two operators
shared_ptr<Vector_Operator> operator-(shared_ptr<Vector_Operator> const op1,
                                      shared_ptr<Vector_Operator> const op2)
{
    return make_shared<Vector_Operator_Difference>(op1, op2);
}

// Product of two operators
shared_ptr<Vector_Operator> operator*(shared_ptr<Vector_Operator> const op1,
                                      shared_ptr<Vector_Operator> const op2)
{
    return make_shared<Vector_Operator_Product>(op1, op2);
}
