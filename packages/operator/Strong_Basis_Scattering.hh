#ifndef Strong_Basis_Scattering_hh
#define Strong_Basis_Scattering_hh

#include <memory>
#include <vector>

#include "Full_Scattering_Operator.hh"

/*
  Applies scattering to a moment representation of the flux
*/
class Strong_Basis_Scattering : public Full_Scattering_Operator
{
public:

    // Constructor
    Strong_Basis_Scattering(std::shared_ptr<Weak_Spatial_Discretization> spatial_discretization,
                            std::shared_ptr<Angular_Discretization> angular_discretization,
                            std::shared_ptr<Energy_Discretization> energy_discretization,
                            Options options = Options());
    
    virtual void check_class_invariants() const override;
    
    virtual std::string description() const override
    {
        return "Strong_Basis_Scattering";
    }
    
private: 

    // Apply within-group and out-of-group scattering
    virtual void apply_full(std::vector<double> &x) const override;
    
    // Apply only within-group scattering
    virtual void apply_coherent(std::vector<double> &x) const override;
};

#endif
