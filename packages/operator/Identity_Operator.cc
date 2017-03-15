#include "Identity_Operator.hh"

using std::vector;

Identity_Operator::
Identity_Operator(int size):
    Vector_Operator(size,
                    size)
{
}

void Identity_Operator::
apply(vector<double> &x) const
{

}

void Identity_Operator::
check_class_invariants() const
{
}
