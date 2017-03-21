#ifndef Sweep_Operator_hh
#define Sweep_Operator_hh

#include <memory>

#include "Square_Vector_Operator.hh"

class Angular_Discretization;
class Energy_Discretization;
class Spatial_Discretization;
class Transport_Discretization;
class XML_Node;

/*
  Transport inverse operator
  
  \Omega \cdot \nabla \psi + \Sigma_t \psi = q
*/
class Sweep_Operator : public Square_Vector_Operator
{
public:

    enum class Sweep_Type
    {
        MOMENT,
        ORDINATE
    };
    
    // Constructor
    Sweep_Operator(Sweep_Type sweep_type,
                   std::shared_ptr<Transport_Discretization> transport_discretization);

    virtual int size() const override
    {
        return size_;
    }
    
    // Include boundary source in sweep
    virtual bool include_boundary_source() const
    {
        return include_boundary_source_;
    }
    virtual void set_include_boundary_source(bool include_source)
    {
        include_boundary_source_ = include_source;
    }

    // Sweep type
    virtual Sweep_Type sweep_type() const
    {
        return sweep_type_;
    }
    
    // Data
    virtual std::shared_ptr<Spatial_Discretization> spatial_discretization() const = 0;
    virtual std::shared_ptr<Angular_Discretization> angular_discretization() const = 0;
    virtual std::shared_ptr<Energy_Discretization> energy_discretization() const = 0;
    virtual std::shared_ptr<Transport_Discretization> transport_discretization() const
    {
        return transport_discretization_;
    }

    virtual void output(XML_Node output_node) const = 0;
    virtual void check_class_invariants() const override = 0;
    
protected:
    
    bool include_boundary_source_;
    Sweep_Type sweep_type_;
    std::shared_ptr<Transport_Discretization> transport_discretization_;
    
private:

    int size_;
    virtual void apply(std::vector<double> &x) const override = 0;
};

#endif
