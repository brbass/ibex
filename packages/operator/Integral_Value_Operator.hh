#ifndef Integral_Value_Operator_hh
#define Integral_Value_Operator_hh

#include "Integration_Mesh.hh"
#include "Vector_Operator.hh"

#include <memory>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;
class Weak_Spatial_Discretization;

class Integral_Value_Operator : public Vector_Operator
{
public:

    Integral_Value_Operator(std::shared_ptr<Integration_Mesh_Options> options,
                            std::shared_ptr<Weak_Spatial_Discretization> spatial,
                            std::shared_ptr<Angular_Discretization> angular,
                            std::shared_ptr<Energy_Discretization> energy);

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
        return "Integral_Value_Operator";
    }

private:

    virtual void apply(std::vector<double> &x) const override;
    
    void get_flux(std::shared_ptr<Integration_Mesh::Cell> const cell,
                  std::vector<double> const &b_val,
                  std::vector<double> const &coeff,
                  std::vector<double> &flux) const;
    
    int row_size_;
    int column_size_;
    
    std::shared_ptr<Integration_Mesh_Options> integration_options_;
    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
    
    mutable bool initialized_;
    mutable std::shared_ptr<Integration_Mesh> mesh_;
};

#endif
   
