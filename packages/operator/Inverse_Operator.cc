#include "Inverse_Operator.hh"

using std::shared_ptr;
using std::string;
using std::vector;

Inverse_Operator::
Inverse_Operator(std::shared_ptr<Vector_Operator> vector_operator):
    Square_Vector_Operator(),
    vector_operator_(vector_operator)
{
    Assert(vector_operator_->square());
    size_ = vector_operator_->column_size();
}

