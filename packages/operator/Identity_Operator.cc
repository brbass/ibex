#include "Identity_Operator.hh"

using std::vector;

Identity_Operator::
Identity_Operator(int size):
    Square_Vector_Operator(),
    size_(size)
{
    check_class_invariants();
}

void Identity_Operator::
apply(vector<double> &x) const
{

}

void Identity_Operator::
check_class_invariants() const
{
    Assert(size_ >= 0);
}
