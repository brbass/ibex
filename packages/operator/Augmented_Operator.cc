#include "Augmented_Operator.hh"

using std::shared_ptr;
using std::vector;

Augmented_Operator::
Augmented_Operator(int number_of_augments,
                   shared_ptr<Vector_Operator> vector_operator,
                   bool zero_out_augments):
    Vector_Operator(),
    zero_out_augments_(zero_out_augments),
    number_of_augments_(number_of_augments),
    vector_operator_(vector_operator)
{
    check_class_invariants();
}

void Augmented_Operator::
apply(vector<double> &x) const
{
    int operator_row_size = vector_operator_->row_size();
    int operator_column_size = vector_operator_->column_size();
    
    vector<double> y(x.begin() + operator_column_size, x.end());
    x.resize(operator_column_size);
    (*vector_operator_)(x);
    
    if (zero_out_augments_)
    {
        x.resize(row_size(), 0);
    }
    else
    {
        x.insert(x.end(), y.begin(), y.end());
    }
}

void Augmented_Operator::
check_class_invariants() const
{
    Assert(vector_operator_);
    Assert(number_of_augments_ >= 0);
}
