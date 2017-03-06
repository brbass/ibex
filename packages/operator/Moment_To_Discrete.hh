#ifndef Moment_To_Discrete_hh
#define Moment_To_Discrete_hh

#include <memory>
#include <vector>

#include "Vector_Operator.hh"

class Angular_Discretization;
class Energy_Discretization;
class Spatial_Discretization;

/*
  Converts moments of the angular flux to the discrete angular flux
*/
class Moment_To_Discrete: public Vector_Operator
{
public:

    // Constructor
    Moment_To_Discrete(std::shared_ptr<Spatial_Discretization> spatial_discretization,
                       std::shared_ptr<Angular_Discretization> angular_discretization,
                       std::shared_ptr<Energy_Discretization> energy_discretization,
                       bool include_dimensional_moments = false);
    
    virtual void check_class_invariants() const override;

private:

    virtual void apply(std::vector<double> &x) const override;

    bool include_dimensional_moments_;
    std::shared_ptr<Spatial_Discretization> spatial_discretization_;
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;
};

#endif
