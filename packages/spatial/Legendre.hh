#ifndef Legendre_hh
#define Legendre_hh

#include <vector>

/*
  One-dimensional Legendre function
  Only used in Legendre_Function for tensor-product functions
*/
class Legendre
{
public:

    Legendre(int order,
             std::vector<double> limits);

    double value(double r) const;
    double d_value(double r) const;
    double dd_value(double r) const;

private:

    double get_local_position(double r) const;
    double get_global_position(double xi) const;
    
    int order_;
    double length_;
    double invlength_;
    std::vector<double> limits_;
};

#endif
