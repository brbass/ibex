#ifndef Identity_Operator_hh
#define Identity_Operator_hh

#include <vector>

#include "Square_Vector_Operator.hh"

/*
  Return the given vector
*/
class Identity_Operator : public Square_Vector_Operator
{
public:

    // Constructor
    Identity_Operator(int size);
    
    virtual void check_class_invariants() const override;

    virtual int size() const override
    {
        return size_;
    }
    
private:

    virtual void apply(std::vector<double> &x) const override;

    int size_;
};

#endif
