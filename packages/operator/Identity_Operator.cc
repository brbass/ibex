#include "Identity_Operator.hh"

using std::vector;

Identity_Operator::
Identity_Operator(int size):
    Square_Vector_Operator(),
    size_(size)
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
