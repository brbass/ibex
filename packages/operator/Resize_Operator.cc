#include "Resize_Operator.hh"

#include "Check.hh"

Resize_Operator::
Resize_Operator(int new_size,
                int original_size):
    new_size_(new_size),
    original_size_(original_size)
{
    check_class_invariants();
}

void Resize_Operator::
apply(std::vector<double> &x) const
{
    x.resize(new_size_, 0);
}

void Resize_Operator::
check_class_invariants() const
{
    Assert(new_size_ >= 0);
    Assert(original_size_ >= 0);
}
