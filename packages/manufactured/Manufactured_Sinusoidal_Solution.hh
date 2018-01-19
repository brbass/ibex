#ifndef Manufactured_Sinusoidal_Solution_hh
#define Manufactured_Sinusoidal_Solution_hh

#include <memory>
#include <vector>

#include "Manufactured_Solution.hh"

class Angular_Discretization;
class Energy_Discretization;

/*
  Solution of the form
  \phi = \phi_0 (1 + u \cos (2\pi f_x x) \cos (2\pi f_y y) \cos (2\pi f_z z) )
  f_d: frequency
  u: relative amplitude
*/
class Manufactured_Sinusoidal_Solution : public Manufactured_Solution
{
public:
    
    Manufactured_Sinusoidal_Solution(std::shared_ptr<Angular_Discretization> angular,
                                     std::shared_ptr<Energy_Discretization> energy,
                                     double const relative_amplitude,
                                     std::vector<double> const &frequency,
                                     std::vector<double> const &solution);
    
    // Solution to problem for all groups and moments (g, m)
    virtual std::vector<double> get_solution(std::vector<double> const &position) const override;
    
    // Gradient of solution to problem for all groups and moments (d, g, m)
    virtual std::vector<double> get_grad_solution(std::vector<double> const &position) const override;
    
private:
    
    double relative_amplitude_;
    std::vector<double> frequency_;
    std::vector<double> solution_;
};
    
#endif
