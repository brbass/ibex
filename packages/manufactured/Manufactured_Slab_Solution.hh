#ifndef Manufactured_Slab_Solution_hh
#define Manufactured_Slab_Solution_hh

#include <memory>
#include <vector>

#include "Manufactured_Solution.hh"

class Angular_Discretization;
class Energy_Discretization;

/*
  Manufactured solution data necessary to create a solid geometry
*/
class Manufactured_Slab_Solution : public Manufactured_Solution
{
public:
    
    Manufactured_Slab_Solution(std::shared_ptr<Angular_Discretization> angular,
                               std::shared_ptr<Energy_Discretization> energy,
                               std::vector<double> interface_positions,
                               std::vector<std::vector<double> > solution);
    
    // Solution to problem for all groups and moments (g, m)
    virtual std::vector<double> get_solution(std::vector<double> const &position) const override;
    
    // Gradient of solution to problem for all groups and moments (d, g, m)
    virtual std::vector<double> get_grad_solution(std::vector<double> const &position) const override;
    
private:

    int number_of_regions_;
    std::vector<double> interface_positions_;
    std::vector<std::vector<double> > solution_;
    std::vector<double> grad_solution_;
};
    
#endif
