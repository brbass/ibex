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
    Square_Vector_Operator();
    
    // Are the input and output size the same?
    virtual bool square() const override
    {
        return true;
    }

    // Input and output size should be the same
    virtual int size() const = 0;
    
    // Output size
    virtual int row_size() const override
    {
        return size();
    }
    
    // Input size
    virtual int column_size() const override
    {
        return size();
    }
    
    virtual void check_class_invariants() const override = 0;

    virtual std::string description() const override = 0;
    
private:
    
    virtual void apply(std::vector<double> &x) const override = 0;
};

#endif
