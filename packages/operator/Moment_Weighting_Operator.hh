#ifndef Moment_Weighting_Operator_hh
#define Moment_Weighting_Operator_hh

#include "Square_Vector_Operator.hh"

#include <memory>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;
class Weak_Spatial_Discretization;

class Moment_Weighting_Operator : public Square_Vector_Operator
{
public:

    Moment_Weighting_Operator(std::shared_ptr<Weak_Spatial_Discretization> spatial,
                              std::shared_ptr<Angular_Discretization> angular,
                              std::shared_ptr<Energy_Discretization> energy);
    
    virtual void check_class_invariants() const override;
    virtual void apply(std::vector<double> &x) const override;

private:

    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
};

#endif
