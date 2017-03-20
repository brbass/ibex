#ifndef Dimensional_Moment_Copy_hh
#define Dimensional_Moment_Copy_hh

#include "Vector_Operator.hh"

#include <memory>

class Angular_Discretization;
class Energy_Discretization;
class Weak_Spatial_Discretization;

class Dimensional_Moment_Copy : public Vector_Operator
{
public:

    // Constructor
    Dimensional_Moment_Copy(std::shared_ptr<Weak_Spatial_Discretization> spatial,
                            std::shared_ptr<Angular_Discretization> angular,
                            std::shared_ptr<Energy_Discretization> energy);
    
    virtual void check_class_invariants() const override;
    
private:

    virtual void apply(std::vector<double> &x) const override;

    // Data
    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
};

#endif