#ifndef Boundary_Source_Toggle_hh
#define Boundary_Source_Toggle_hh

#include <memory>

#include "Sweep_Operator.hh"

/*
  Transport inverse operator
  
  \Omega \cdot \nabla \psi + \Sigma_t \psi = q
*/
class Boundary_Source_Toggle : public Sweep_Operator
{
public:

    // Constructor
    Boundary_Source_Toggle(bool local_include_boundary_source,
                           std::shared_ptr<Sweep_Operator> sweep);
    
    // Data
    virtual std::shared_ptr<Spatial_Discretization> spatial_discretization() const override
    {
        return sweep_->spatial_discretization();
    }
    virtual std::shared_ptr<Angular_Discretization> angular_discretization() const override
    {
        return sweep_->angular_discretization();
    }
    virtual std::shared_ptr<Energy_Discretization> energy_discretization() const override
    {
        return sweep_->energy_discretization();
    }
    virtual std::shared_ptr<Transport_Discretization> transport_discretization() const override
    {
        return sweep_->transport_discretization();
    }
    
    virtual void output(XML_Node output_node) const override;
    virtual void check_class_invariants() const override
    {
        sweep_->check_class_invariants();
    }
    virtual std::string description() const override;
    
private:

    // Data
    bool local_include_boundary_source_;
    std::shared_ptr<Sweep_Operator> sweep_;
    
    virtual void apply(std::vector<double> &x) const override;
};

#endif
