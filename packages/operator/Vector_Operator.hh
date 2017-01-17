#ifndef Vector_Operator_hh
#define Vector_Operator_hh

#include <vector>

#include "Check.hh"

/*
  Pure virtual class to represent a vector operator
*/
class Vector_Operator
{
public:

    // Constructor
    Vector_Operator(int row_size,
                    int column_size);

    // Apply the operator
    std::vector<double> &operator()(std::vector<double> &x)
    {
        Check(x.size() == column_size_);
        
        apply(x);
        
        Check(x.size() == row_size_);
        
        return x;
    }

    // Output size
    virtual int row_size() const
    {
        return row_size_;
    }

    // Input size
    virtual int column_size() const
    {
        return column_size_;
    }

    // Are the input and output size the same?
    virtual bool square() const
    {
        return (row_size_ == column_size_);
    }

    virtual void check_class_invariants() const = 0;
    
private:
    
    virtual void apply(std::vector<double> &x) const = 0;
    
    int row_size_;
    int column_size_;
};

#endif
