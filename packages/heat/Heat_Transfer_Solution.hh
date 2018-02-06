#ifndef Heat_Transfer_Solution_hh
#define Heat_Transfer_Solution_hh

#include <memory>
#include <vector>

class Weak_Spatial_Discretization;

class Heat_Transfer_Solution
{
public:
    
    Heat_Transfer_Solution(std::shared_ptr<Weak_Spatial_Discretization> spatial,
                           std::vector<double> const &coefficients);

    double solution(std::vector<double> const &position) const;
    
private:

    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
    std::vector<double> coefficients_;
};

#endif
