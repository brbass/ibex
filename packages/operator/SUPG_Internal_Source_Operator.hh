#ifndef SUPG_Internal_Source_Operator_hh
#define SUPG_Internal_Source_Operator_hh

#include <memory>
#include <vector>

#include "Square_Vector_Operator.hh"

class Angular_Discretization;
class Energy_Discretization;
class Spatial_Discretization;

/*
  Return the internal source
*/
class SUPG_Internal_Source_Operator : public Vector_Operator
{
public:
    
    // Constructor
    SUPG_Internal_Source_Operator(std::shared_ptr<Spatial_Discretization> spatial_discretization,
                                  std::shared_ptr<Angular_Discretization> angular_discretization,
                                  std::shared_ptr<Energy_Discretization> energy_discretization);
    
    virtual int row_size() const override
    {
        return row_size_;
    }
    virtual int column_size() const override
    {
        return column_size_;
    }
    
    virtual void check_class_invariants() const override;

    virtual std::string description() const override
    {
        return "SUPG_Internal_Source_Operator";
    }
    
private:
    
    virtual void apply(std::vector<double> &x) const override;

    // Data
    int row_size_;
    int column_size_;
    std::shared_ptr<Spatial_Discretization> spatial_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
    
};

#endif
