#ifndef Sweep_Operator_hh
#define Sweep_Operator_hh

#include <memory>

#include "Vector_Operator.hh"

class Angular_Discretization;
class Energy_Discretization;
class Spatial_Discretization;
class Transport_Discretization;

#include "pugixml.hh"

using std::shared_ptr;

/*
  Transport inverse operator
  
  \Omega \cdot \nabla \psi + \Sigma_t \psi = q
*/
class Sweep_Operator : public Vector_Operator
{
public:

    enum class Sweep_Type
    {
        MOMENT,
        ORDINATE
    };
    
    // Constructor
    Sweep_Operator(Sweep_Type sweep_type,
                   shared_ptr<Transport_Discretization> transport_discretization);
    
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
    virtual shared_ptr<Spatial_Discretization> spatial_discretization() const = 0;
    virtual shared_ptr<Angular_Discretization> angular_discretization() const = 0;
    virtual shared_ptr<Energy_Discretization> energy_discretization() const = 0;
    virtual shared_ptr<Transport_Discretization> transport_discretization() const
    {
        return transport_discretization_;
    }

    virtual void output(pugi::xml_node output_node) const = 0;
    virtual void check_class_invariants() const override = 0;
    
protected:
    
    bool include_boundary_source_;
    Sweep_Type sweep_type_;
    shared_ptr<Transport_Discretization> transport_discretization_;
    
private:
    
    virtual void apply(vector<double> &x) const override = 0;
};

#endif
