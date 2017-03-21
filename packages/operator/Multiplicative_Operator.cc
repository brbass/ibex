#include "Multiplicative_Operator.hh"

using std::vector;

Multiplicative_Operator::
Multiplicative_Operator(int size,
                        double scalar):
    Square_Vector_Operator(),
    size_(size),
    scalar_(scalar)
{
    check_class_invariants();
}

void Multiplicative_Operator::
apply(vector<double> &x) const
{
    for (int i = 0; i < size_; ++i)
    {
        x[i] *= scalar_;
    }
}

void Multiplicative_Operator::
check_class_invariants() const
{
}
