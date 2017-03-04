#ifndef Scattering_hh
#define Scattering_hh

#include <memory>
#include <vector>

#include "Scattering_Operator.hh"

/*
  Applies scattering to a moment representation of the flux
*/
class Scattering : public Scattering_Operator
{
public:

    // Constructor
    Scattering(std::shared_ptr<Spatial_Discretization> spatial_discretization,
               std::shared_ptr<Angular_Discretization> angular_discretization,
               std::shared_ptr<Energy_Discretization> energy_discretization,
               Scattering_Type scattering_type = Scattering_Type::FULL);
    
    virtual void check_class_invariants() const override;
    
private: 

    // Apply within-group and out-of-group scattering
    virtual void apply_full(std::vector<double> &x) const override;
    
    // Apply only within-group scattering
    virtual void apply_coherent(std::vector<double> &x) const override;
};

#endif
