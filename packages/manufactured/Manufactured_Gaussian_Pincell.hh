#ifndef Manufactured_Gaussian_Pincell_hh
#define Manufactured_Gaussian_Pincell_hh

#include <memory>
#include <vector>

#include "Manufactured_Solution.hh"

class Angular_Discretization;
class Distance;
class Energy_Discretization;

/*
  Manufactured solution of the form
  phi_0 exp(a sqrt(x^2 + y^2 + z^2)) (1 + b * exp(c * (sqrt(x^2 + y^2 + z^2) + d)^2))
  The constants a, b, c, d are dependent on group and moment
*/
class Manufactured_Gaussian_Pincell : public Manufactured_Solution
{
public:
    
    Manufactured_Gaussian_Pincell(std::shared_ptr<Angular_Discretization> angular,
                                  std::shared_ptr<Energy_Discretization> energy,
                                  std::vector<double> solution,
                                  std::vector<double> aval,
                                  std::vector<double> bval,
                                  std::vector<double> cval,
                                  std::vector<double> dval);
    
    // Solution to problem for all groups and moments (g, m)
    virtual std::vector<double> get_solution(std::vector<double> const &position) const override;
    
    // Gradient of solution to problem for all groups and moments (d, g, m)
    virtual std::vector<double> get_grad_solution(std::vector<double> const &position) const override;
    
private:

    std::vector<double> solution_;
    std::vector<double> aval_;
    std::vector<double> bval_;
    std::vector<double> cval_;
    std::vector<double> dval_;
    std::vector<double> origin_;
    std::shared_ptr<Distance> distance_;
};
    
#endif
