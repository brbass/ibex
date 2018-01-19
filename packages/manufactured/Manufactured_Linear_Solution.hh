#ifndef Manufactured_Linear_Solution_hh
#define Manufactured_Linear_Solution_hh

#include <memory>
#include <vector>

#include "Manufactured_Solution.hh"

class Angular_Discretization;
class Energy_Discretization;

/*
  Manufactured solution of the form
  phi_0 (1 + \sum_d s_d (x_d -x_{d,0}))
*/
class Manufactured_Linear_Solution : public Manufactured_Solution
{
public:
    
    Manufactured_Linear_Solution(std::shared_ptr<Angular_Discretization> angular,
                                 std::shared_ptr<Energy_Discretization> energy,
                                 std::vector<double> origin,
                                 std::vector<double> slope,
                                 std::vector<double> solution);
    
    // Solution to problem for all groups and moments (g, m)
    virtual std::vector<double> get_solution(std::vector<double> const &position) const override;
    
    // Gradient of solution to problem for all groups and moments (d, g, m)
    virtual std::vector<double> get_grad_solution(std::vector<double> const &position) const override;
    
private:

    std::vector<double> origin_;
    std::vector<double> slope_;
    std::vector<double> solution_;
};
    
#endif
