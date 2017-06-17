#ifndef Combined_SUPG_Fission_hh
#define Combined_SUPG_Fission_hh

#include <memory>
#include <vector>

#include "Combined_SUPG_Operator.hh"

/*
  Applies scattering to a moment representation of the flux
*/
class Combined_SUPG_Fission : public Combined_SUPG_Operator
{
public:

    // Constructor
    Combined_SUPG_Fission(std::shared_ptr<Spatial_Discretization> spatial_discretization,
                          std::shared_ptr<Angular_Discretization> angular_discretization,
                          std::shared_ptr<Energy_Discretization> energy_discretization,
                          Options options = Options());
    
    virtual void check_class_invariants() const override;
    
    virtual std::string description() const override
    {
        return "Combined_SUPG_Fission";
    }
    
private: 

    // Apply within-group and out-of-group scattering
    virtual void apply_full(std::vector<double> &x) const override;
    
    // Apply only within-group scattering
    virtual void apply_coherent(std::vector<double> &x) const override;
};

#endif
