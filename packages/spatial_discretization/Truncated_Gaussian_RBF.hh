#ifndef Truncated_Gaussian_RBF_hh
#define Truncated_Gaussian_RBF_hh

#include "Local_RBF.hh"

class Truncated_Gaussian_RBF : public Local_RBF
{
public:
    
    // Constructor
    Truncated_Gaussian_RBF(double radius = 5);

    // Distance from center the function is nonzero
    virtual double radius() const
    {
        return radius_;
    }
    
    // Value of basis function
    virtual double value(double r) const override;
    
    // Derivative of basis function
    virtual double d_value(double r) const override;
    
    // Second derivative of the basis function
    virtual double dd_value(double r) const override;

    // Description of RBF
    virtual std::string description() const override
    {
        return "truncated_gaussian";
    }

private:
    
    double radius_;
};

#endif
