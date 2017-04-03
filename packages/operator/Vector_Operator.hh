#ifndef Vector_Operator_hh
#define Vector_Operator_hh

#include <string>
#include <vector>

#include "Check.hh"

/*
  Pure virtual class to represent a vector operator
*/
class Vector_Operator
{
public:

    // Constructor
    Vector_Operator();

    // Apply the operator
    std::vector<double> &operator()(std::vector<double> &x)
    {
        Check(x.size() == column_size());
        
        apply(x);
        
        Check(x.size() == row_size());
        
        return x;
    }

    // Output size
    virtual int row_size() const = 0;

    // Input size
    virtual int column_size() const = 0;
    
    // Are the input and output size the same?
    virtual bool square() const
    {
        return (row_size() == column_size());
    }
    
    virtual void check_class_invariants() const = 0;
    virtual std::string description() const = 0;
    
private:
    
    virtual void apply(std::vector<double> &x) const = 0;
};

#endif
