#ifndef Null_Operator_hh
#define Null_Operator_hh

#include <vector>

#include "Vector_Operator.hh"

using std::vector;

/*
  Zeroes out any vector given it
*/
class Null_Operator : public Vector_Operator
{
public:

    // Constructor
    Null_Operator(int size);

    virtual void check_class_invariants() const override
    {
    }
    
private:
    
    virtual void apply(vector<double> &x) const override;
};

#endif
