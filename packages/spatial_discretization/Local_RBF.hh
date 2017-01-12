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
    virtual double radius() const = 0;
    
    // Value of basis function
    virtual double basis(double r) const = 0;
    
    // Derivative of basis function
    virtual double d_basis(double r) const = 0;
    
    // Second derivative of the basis function
    virtual double dd_basis(double r) const = 0;

    // Description of RBF
    virtual std::string description() const = 0;
};

#endif
