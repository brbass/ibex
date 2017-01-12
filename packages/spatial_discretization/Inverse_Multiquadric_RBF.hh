#ifndef Inverse_Multiquadric_RBF_hh
#define Inverse_Multiquadric_RBF_hh

#include "Global_RBF.hh"

/*
  Radial basis function of the form exp(-c^2 r^2)
*/
class Inverse_Multiquadric_RBF : public Global_RBF
{
public:

    // Constructor
    Inverse_Multiquadric_RBF();

    // Value of basis function
    virtual double basis(double r) const override;
    
    // Derivative of basis function
    virtual double d_basis(double r) const override;
    
    // Second derivative of the basis function
    virtual double dd_basis(double r) const override;

    virtual std::string description() const override
    {
        return "inverse_multiquadric";
    }
};

#endif
