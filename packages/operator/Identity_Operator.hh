#ifndef Identity_Operator_hh
#define Identity_Operator_hh

#include <vector>

#include "Vector_Operator.hh"

/*
  Return the given vector
*/
class Identity_Operator : public Vector_Operator
{
public:

    // Constructor
    Identity_Operator(int size);
    
    virtual void check_class_invariants() const override;
    
private:

    virtual void apply(std::vector<double> &x) const override;
};

#endif
