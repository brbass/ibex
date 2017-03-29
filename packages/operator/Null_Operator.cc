#include "Null_Operator.hh"

using std::vector;

Null_Operator::
Null_Operator(int size):
    Square_Vector_Operator(),
    size_(size)
{
    check_class_invariants();
}

void Null_Operator::
apply(vector<double> &x) const
{
    x.assign(x.size(), 0);
}

void Null_Operator::
check_class_invariants() const
{
}
