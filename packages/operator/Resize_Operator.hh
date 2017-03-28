#ifndef Resize_Operator_hh
#define Resize_Operator_hh

#include <memory>
#include <vector>

#include "Vector_Operator.hh"

/*
  Resize vector operator and put zeros in if size is larger
*/
class Resize_Operator: public Vector_Operator
{
public:

    // Constructor
    Resize_Operator(int new_size,
                    int original_size);
    
    virtual void check_class_invariants() const override;
    
    virtual int row_size() const override
    {
        return new_size_;
    }
    virtual int column_size() const override
    {
        return original_size_;
    }
    
private:

    virtual void apply(std::vector<double> &x) const override;
    
    int new_size_;
    int original_size_;
};

#endif
