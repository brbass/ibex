#ifndef RBF_hh
#define RBF_hh

#include <string>

/*
  Pure virtual class to represent a radial basis function
  
  r/s: nondimensional distance
  ds: derivative of nondimensional distance
  dds: second derivative of nondimensional distance
*/
class RBF
{
public:

    enum class Range
    {
        LOCAL,
        GLOBAL
    };
    
    // Constructor
    RBF();

    // Range of function
    virtual Range range() const = 0;
    
    // Distance from center the function is nonzero
    virtual double radius() const = 0;
    
    // Value of basis function
    virtual double value(double r) const = 0;
    
    // Derivative of basis function
    virtual double d_value(double r) const = 0;
    virtual double d_value(double s,
                           double ds) const;
    
    // Second derivative of the basis function
    virtual double dd_value(double r) const = 0;
    virtual double dd_value(double s,
                            double ds,
                            double dds) const;
    
    virtual std::string description() const = 0;
};

#endif
