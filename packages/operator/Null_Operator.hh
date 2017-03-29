#ifndef Null_Operator_hh
#define Null_Operator_hh

#include <vector>

#include "Square_Vector_Operator.hh"

/*
  Zeroes out any vector given it
*/
class Null_Operator : public Square_Vector_Operator
{
public:

    // Constructor
    Null_Operator(int size);

    virtual int size() const override
    {
        return size_;
    }
    
    virtual void check_class_invariants() const override;
    
private:

    int size_;
    
    virtual void apply(std::vector<double> &x) const override;
};

#endif
