#ifndef Local_RBF_hh
#define Local_RBF_hh

#include "RBF.hh"

class Local_RBF : public RBF
{
public:
    // Constructor
    Local_RBF();

    // Range of function
    virtual Range range() const override
    {
        return Range::LOCAL;
    }
    
    // Distance from center the function is nonzero
    virtual double radius() const override = 0;
    
    // Value of basis function
    virtual double value(double r) const override = 0;
    
    // Derivative of basis function
    virtual double d_value(double r) const override = 0;
    
    // Second derivative of the basis function
    virtual double dd_value(double r) const override = 0;

    // Description of RBF
    virtual std::string description() const override = 0;
};

#endif
