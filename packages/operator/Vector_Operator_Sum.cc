#include "Vector_Operator_Sum.hh"

Vector_Operator_Sum::
Vector_Operator_Sum(shared_ptr<Vector_Operator> op1,
                    shared_ptr<Vector_Operator> op2):
    Vector_Operator(op1->row_size(),
                    op1->column_size()),
    op1_(op1),
    op2_(op2)
{
}

void Vector_Operator_Sum::
apply(vector<double> &x) const
{
    vector<double> y(x);

    (*op1_)(x);
    (*op2_)(y);

    for (int i = 0; i < row_size(); ++i)
    {
        x[i] += y[i];
    }
}

void Vector_Operator_Sum::
check_class_invariants() const
{
    Assert(op1_);
    Assert(op2_);
}
