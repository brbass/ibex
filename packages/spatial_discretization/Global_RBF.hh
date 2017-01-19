#ifndef Global_RBF_hh
#define Global_RBF_hh

#include "RBF.hh"

class Global_RBF : public RBF
{
public:

    // Constructor
    Global_RBF();

    // Range of function
    virtual Range range() const override
    {
        return Range::GLOBAL;
    }
    
    // Distance from center the function is nonzero
    virtual double radius() const override;
    
    // Value of basis function
    virtual double value(double r) const = 0;
    
    // Derivative of basis function
    virtual double d_value(double r) const = 0;
    
    // Second derivative of the basis function
    virtual double dd_value(double r) const = 0;

    // Description of RBF
    virtual std::string description() const = 0;
};

#endif
