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
                       std::shared_ptr<Energy_Discretization> energy_discretization,
                       bool include_dimensional_moments = false);
    
    virtual void check_class_invariants() const override;
    
    virtual int row_size() const override
    {
        return row_size_;
    }
    virtual int column_size() const override
    {
        return column_size_;
    }
    
private:
    
    virtual void apply(std::vector<double> &x) const override;

    bool include_dimensional_moments_;
    int row_size_;
    int column_size_;
    std::shared_ptr<Spatial_Discretization> spatial_discretization_;
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;
};

#endif
