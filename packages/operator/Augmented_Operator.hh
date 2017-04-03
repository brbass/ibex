#ifndef Augmented_Operator_hh
#define Augmented_Operator_hh

#include <memory>

#include "Vector_Operator.hh"

/*
  Wraps a Vector_Operator to allow unmodified augments at the end of the input vector
*/
class Augmented_Operator : public Vector_Operator
{
public:

    // Constructor
    Augmented_Operator(int number_of_augments,
                       std::shared_ptr<Vector_Operator> vector_operator,
                       bool zero_out_augments = false);
    
    virtual void check_class_invariants() const override;

    virtual int row_size() const override
    {
        return vector_operator_->row_size() + number_of_augments_;
    }
    virtual int column_size() const override
    {
        return vector_operator_->column_size() + number_of_augments_;
    }
    virtual std::string description() const override;
    
private:
    
    virtual void apply(std::vector<double> &x) const override;

    bool zero_out_augments_;
    int number_of_augments_;
    std::shared_ptr<Vector_Operator> vector_operator_;
};

#endif
