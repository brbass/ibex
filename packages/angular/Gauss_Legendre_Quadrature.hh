#ifndef Gauss_Legendre_Quadrature_hh
#define Gauss_Legendre_Quadrature_hh

#include <vector>

#include "Angular_Discretization.hh"

/*
  Holds Gauss Legendre quadrature
*/
class Gauss_Legendre_Quadrature : public Angular_Discretization
{
public:

    // Constructor
    Gauss_Legendre_Quadrature(int dimension,
                              int number_of_moments,
                              int number_of_ordinates);

    // Return direction
    virtual std::vector<double> const &direction(int ord) const override
    {
        return directions_[ord];
    }
    
    // Return Gauss-Legendre ordinates
    virtual std::vector<double> const &ordinates() const override
    {
        return ordinates_;
    }
    
    // Return Gauss-Legendre weights
    virtual std::vector<double> const &weights() const override
    {
        return weights_;
    }

    // Check class invariants
    virtual void check_class_invariants() const override;

    // Output data to XML file
    virtual void output(XML_Node output_node) const override;

    virtual int reflect_ordinate(int o,
                                 std::vector<double> const &/*normal*/) const override
    {
        return number_of_ordinates_ - o - 1;
    }
    
private:

    std::vector<double> ordinates_;
    std::vector<double> weights_;
    std::vector<std::vector<double> > directions_;
};

#endif
