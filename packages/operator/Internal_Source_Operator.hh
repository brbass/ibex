#ifndef Internal_Source_Operator_hh
#define Internal_Source_Operator_hh

#include <memory>
#include <vector>

#include "Square_Vector_Operator.hh"

class Angular_Discretization;
class Energy_Discretization;
class Spatial_Discretization;

/*
  Return the internal source
*/
class Internal_Source_Operator : public Square_Vector_Operator
{
public:
    
    // Constructor
    Internal_Source_Operator(std::shared_ptr<Spatial_Discretization> spatial_discretization,
            std::shared_ptr<Angular_Discretization> angular_discretization,
            std::shared_ptr<Energy_Discretization> energy_discretization);
    
    // Input and output size should be the same
    virtual int size() const override
    {
        return size_;
    }
    
    virtual void check_class_invariants() const override;
    
private:
    
    virtual void apply(std::vector<double> &x) const override;

    // Data
    int size_;
    std::shared_ptr<Spatial_Discretization> spatial_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
    
};

#endif
