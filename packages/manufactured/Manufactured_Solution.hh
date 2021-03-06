#ifndef Manufactured_Solution_hh
#define Manufactured_Solution_hh

#include <memory>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;

/*
  Manufactured solution data necessary to create a solid geometry
*/
class Manufactured_Solution
{
public:
    
    Manufactured_Solution(std::shared_ptr<Angular_Discretization> angular,
                          std::shared_ptr<Energy_Discretization> energy);
    
    // Solution to problem for all groups and moments (g, m)
    virtual std::vector<double> get_solution(std::vector<double> const &position) const = 0;
    
    // Gradient of solution to problem for all groups and moments (d, g, m)
    virtual std::vector<double> get_grad_solution(std::vector<double> const &position) const = 0;
    
    // Get source
    virtual std::vector<double> get_source(std::vector<double> const &solution,
                                           std::vector<double> const &grad_solution,
                                           std::vector<double> const &sigma_t,
                                           std::vector<double> const &sigma_s) const;
    
protected:

    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;

private:
    
    std::vector<int> streaming_size_;
    std::vector<std::vector<int> > streaming_indices_;
    std::vector<std::vector<double> > streaming_coefficients_;
    std::vector<int> expected_zero_moments_;
};
    
#endif
