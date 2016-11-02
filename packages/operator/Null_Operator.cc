#include "Null_Operator.hh"

Null_Operator::
Null_Operator(int size):
    Vector_Operator(size,
                    size)
{
}

void Null_Operator::
apply(vector<double> &x) const
{
    x.assign(x.size(), 0);
}
