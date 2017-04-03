#ifndef Multiplicative_Operator_hh
#define Multiplicative_Operator_hh

#include <vector>

#include "Square_Vector_Operator.hh"

/*
  Return the given vector
*/
class Multiplicative_Operator : public Square_Vector_Operator
{
public:

    // Constructor
    Multiplicative_Operator(int size,
                            double scalar);
    
    void set_scalar(double scalar)
    {
        scalar_ = scalar;
    }

    virtual int size() const override
    {
        return size_;
    }
    
    virtual void check_class_invariants() const override;
    
    virtual std::string description() const override
    {
        return "Multiplicative_Operator";
    }
    
private:

    int size_;
    double scalar_;
    
    virtual void apply(std::vector<double> &x) const override;
};

#endif
