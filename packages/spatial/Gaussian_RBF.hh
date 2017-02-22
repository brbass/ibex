#ifndef Gaussian_RBF_hh
#define Gaussian_RBF_hh

#include "Global_RBF.hh"

/*
  Radial basis function of the form exp(-c^2 r^2)
*/
class Gaussian_RBF : public Global_RBF
{
public:

    // Constructor
    Gaussian_RBF();
    
    // Value of basis function
    virtual double value(double r) const override;
    
    // Derivative of basis function
    virtual double d_value(double r) const override;
    
    // Second derivative of the basis function
    virtual double dd_value(double r) const override;

    virtual std::string description() const override
    {
        return "gaussian";
    }
};

#endif
