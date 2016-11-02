#ifndef Augmented_Operator_hh
#define Augmented_Operator_hh

#include <memory>

#include "Vector_Operator.hh"

using std::shared_ptr;

/*
  Wraps a Vector_Operator to allow unmodified augments at the end of the input vector
*/
class Augmented_Operator : public Vector_Operator
{
public:

    // Constructor
    Augmented_Operator(unsigned number_of_augments,
                       shared_ptr<Vector_Operator> vector_operator,
                       bool zero_out_augments = false);
    
    virtual void check_class_invariants() const override;

private:
    
    virtual void apply(vector<double> &x) const override;

    bool zero_out_augments_;
    int number_of_augments_;
    shared_ptr<Vector_Operator> vector_operator_;
};

#endif
