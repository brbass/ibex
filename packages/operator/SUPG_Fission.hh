#ifndef SUPG_Fission_hh
#define SUPG_Fission_hh

#include <memory>
#include <vector>

#include "Scattering_Operator.hh"

/*
  Applies fission to a moment representation of the flux
*/
class SUPG_Fission : public Scattering_Operator
{
public:

    // Constructor
    SUPG_Fission(std::shared_ptr<Spatial_Discretization> spatial_discretization,
                 std::shared_ptr<Angular_Discretization> angular_discretization,
                 std::shared_ptr<Energy_Discretization> energy_discretization,
                 Options options = Options());
    
    virtual void check_class_invariants() const override;
    
private: 
    
    // Apply within-group and out-of-group fission
    virtual void apply_full(std::vector<double> &x) const override;
    
    // Apply only within-group fission
    virtual void apply_coherent(std::vector<double> &x) const override;
    
    // Regular fission with nu, chi and sigma_f
    void group_full(std::vector<double> &x) const;
    void group_coherent(std::vector<double> &x) const;
    
    // Scattering-like fission (group to group) with all info in sigma_f
    void group_to_group_full(std::vector<double> &x) const;
    void group_to_group_coherent(std::vector<double> &x) const;
};

#endif
