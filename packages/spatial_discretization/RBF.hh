#ifndef RBF_hh
#define RBF_hh

#include <string>

using std::string;

/*
  Pure virtual class to represent a radial basis function
  
  r/s: nondimensional distance
  ds: derivative of nondimensional distance
  dds: second derivative of nondimensional distance
*/
class RBF
{
public:

    // Constructor
    RBF();

    // Value of basis function
    virtual double basis(double r) const = 0;
    
    // Derivative of basis function
    virtual double d_basis(double r) const = 0;
    virtual double d_basis(double s,
                           double ds) const;
    
    // Second derivative of the basis function
    virtual double dd_basis(double r) const = 0;
    virtual double dd_basis(double s,
                            double ds,
                            double dds) const;

    virtual string description() const = 0;

    virtual bool derivative_available(int derivative) const
    {
        return derivative < 3 ? true : false;
    }
};

#endif
