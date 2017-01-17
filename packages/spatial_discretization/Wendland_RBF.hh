#ifndef Wendland_RBF_hh
#define Wendland_RBF_hh

#include "Local_RBF.hh"

/*
  Local polynomial radial basis function
*/
class Wendland_RBF : public Local_RBF
{
public:

    // Constructor
    Wendland_RBF(int order);
    
    // Distance from center the function is nonzero
    virtual double radius() const
    {
        return 1.;
    }
    
    // Value of basis function
    virtual double basis(double r) const override;
    
    // Derivative of basis function
    virtual double d_basis(double r) const override;
    
    // Second derivative of the basis function
    virtual double dd_basis(double r) const override;

    virtual std::string description() const override
    {
        return "wendland" + std::to_string(order_);
    }
    
private:

    int order_;
};

#endif
