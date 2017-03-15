#include "Vector_Operator_Difference.hh"

using std::shared_ptr;
using std::vector;

Vector_Operator_Difference::
Vector_Operator_Difference(shared_ptr<Vector_Operator> op1,
                           shared_ptr<Vector_Operator> op2):
    Vector_Operator(op1->row_size(),
                    op1->column_size()),
    op1_(op1),
    op2_(op2)
{
}

void Vector_Operator_Difference::
apply(vector<double> &x) const
{
    vector<double> y(x);
    
    (*op1_)(x);
    (*op2_)(y);

    for (int i = 0; i < row_size(); ++i)
    {
        x[i] -= y[i];
    }
}

void Vector_Operator_Difference::
check_class_invariants() const
{
    Assert(op1_);
    Assert(op2_);

    // Check that input and output sizes are the same
    int row1 = op1_->row_size();
    int row2 = op2_->row_size();
    int col1 = op1_->column_size();
    int col2 = op2_->column_size();
    Assert(row1 == row2);
    Assert(col1 == col2);
}
