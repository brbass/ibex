#ifndef Discrete_To_Moment_hh
#define Discrete_To_Moment_hh

#include <memory>
#include <vector>

#include "Vector_Operator.hh"

class Angular_Discretization;
class Energy_Discretization;
class Spatial_Discretization;

/*
  Converts the discrete angular flux to moments of the  angular flux
*/
class Discrete_To_Moment: public Vector_Operator
{
public:

    // Constructor
    Discrete_To_Moment(std::shared_ptr<Spatial_Discretization> spatial_discretization,
                       std::shared_ptr<Angular_Discretization> angular_discretization,
                       std::shared_ptr<Energy_Discretization> energy_discretization);
    
    virtual void check_class_invariants() const override;
    
private:
    
    virtual void apply(std::vector<double> &x) const override;

    std::shared_ptr<Spatial_Discretization> spatial_discretization_;
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;
};

#endif
