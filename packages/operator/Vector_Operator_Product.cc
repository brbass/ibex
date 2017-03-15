#include "Vector_Operator_Product.hh"

using std::shared_ptr;
using std::vector;

Vector_Operator_Product::
Vector_Operator_Product(shared_ptr<Vector_Operator> op1,
                        shared_ptr<Vector_Operator> op2):
    Vector_Operator(op1->row_size(),
                    op2->column_size()),
    op1_(op1),
    op2_(op2)
{
}

void Vector_Operator_Product::
apply(vector<double> &x) const
{
    (*op2_)(x);
    (*op1_)(x);
}

void Vector_Operator_Product::
check_class_invariants() const
{
    Assert(op1_);
    Assert(op2_);
    
    // Check that input of op1 is of the size of the output of op2
    int row2 = op2_->row_size();
    int col1 = op1_->column_size();
    Assert(row2 == col1);
}}
