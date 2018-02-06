#ifndef Heat_Transfer_Solve_hh
#define Heat_Transfer_Solve_hh

#include <memory>

class Heat_Transfer_Integration;
class Heat_Transfer_Solution;
class Weak_Spatial_Discretization;

class Heat_Transfer_Solve
{
public:

    Heat_Transfer_Solve(std::shared_ptr<Heat_Transfer_Integration> integration,
                        std::shared_ptr<Weak_Spatial_Discretization> spatial);
    
    std::shared_ptr<Heat_Transfer_Solution> solve();
    
private:

    std::shared_ptr<Heat_Transfer_Integration> integration_;
    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
};

#endif
