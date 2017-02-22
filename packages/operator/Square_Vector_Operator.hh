#ifndef Square_Vector_Operator_hh
#define Square_Vector_Operator_hh

#include <vector>

#include "Vector_Operator.hh"

/*
  Pure virtual class to represent a vector operator
*/
class Square_Vector_Operator : public Vector_Operator
{
public:
    
    // Constructor
    Square_Vector_Operator(int size);
    
    // Are the input and output size the same?
    virtual bool square() const override
    {
        return true;
    }

    virtual void check_class_invariants() const = 0;
    
private:
    
    virtual void apply(std::vector<double> &x) const = 0;
};

#endif
