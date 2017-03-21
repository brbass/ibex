#include "Null_Operator.hh"

using std::vector;

Null_Operator::
Null_Operator(int size):
    Square_Vector_Operator(),
    size_(size)
{
}

void Null_Operator::
apply(vector<double> &x) const
{
    x.assign(x.size(), 0);
}
