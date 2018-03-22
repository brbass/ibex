#ifndef Manufactured_Sin_Data_hh
#define Manufactured_Sin_Data_hh

#include <memory>
#include <vector>

#include "Manufactured_Heat_Data.hh"

/*
  Represents the function t(x,y,z) = t1 + Sin(t2 d(x,y,z)^2),
  where d(x,y,z) is the Cartesian distance.
  The conduction is of the form k(x, y, z) = k1 + k2 * d(x,y,z)^2
  The convection is of the form h(x, y, z) = h1 + h2 * sin(h3 * x * y)
*/
class Manufactured_Sin_Data : public Manufactured_Heat_Data
{
public:
    
    Manufactured_Sin_Data(int dimension,
                          std::vector<std::vector<double> > limits,
                          std::vector<double> sol_coeff,
                          std::vector<double> cond_coeff,
                          std::vector<double> conv_coeff);
    
    // Functions to be implemented for each case
    virtual double conduction(std::vector<double> const &position) const override;
    virtual std::vector<double> gradient_conduction(std::vector<double> const &position) const override;
    virtual double convection(std::vector<double> const &position) const override;
    virtual double solution(std::vector<double> const &position) const override;
    virtual std::vector<double> gradient_solution(std::vector<double> const &position) const override;
    virtual double laplacian_solution(std::vector<double> const &position) const override;
    
private:
    
    std::vector<double> sol_coeff_;
    std::vector<double> cond_coeff_;
    std::vector<double> conv_coeff_;
};

#endif
