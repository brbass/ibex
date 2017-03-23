#ifndef Wendland1_RBF_hh
#define Wendland1_RBF_hh

#include "Local_RBF.hh"

/*
  Local polynomial radial basis function
*/
class Wendland1_RBF : public Local_RBF
{
public:

    // Constructor
    Wendland1_RBF(int order);
    
    // Distance from center the function is nonzero
    virtual double radius() const
    {
        return 1.;
    }
    
    // Value of basis function
    virtual double value(double r) const override;
    
    // Derivative of basis function
    virtual double d_value(double r) const override;
    
    // Second derivative of the basis function
    virtual double dd_value(double r) const override;

    virtual std::string description() const override
    {
        return "wendland" + std::to_string(order_);
    }
    
private:

    int order_;
};

#endif
