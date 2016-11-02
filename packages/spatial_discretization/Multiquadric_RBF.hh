#ifndef Multiquadric_RBF_hh
#define Multiquadric_RBF_hh

#include "RBF.hh"

/*
  Radial basis function of the form exp(-c^2 r^2)
*/
class Multiquadric_RBF : public RBF
{
public:

    // Constructor
    Multiquadric_RBF();

    // Value of basis function
    virtual double basis(double r) const override;
    
    // Derivative of basis function
    virtual double d_basis(double r) const override;
    
    // Second derivative of the basis function
    virtual double dd_basis(double r) const override;

    virtual string description() const override
    {
        return "multiquadric";
    }
};

#endif
