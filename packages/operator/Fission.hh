#ifndef Fission_hh
#define Fission_hh

#include <memory>
#include <vector>

#include "Scattering_Operator.hh"

using std::shared_ptr;
using std::vector;

/*
  Applies fission to a moment representation of the flux
*/
class Fission : public Scattering_Operator
{
public:

    // Constructor
    Fission(shared_ptr<Spatial_Discretization> spatial_discretization,
            shared_ptr<Angular_Discretization> angular_discretization,
            shared_ptr<Energy_Discretization> energy_discretization,
            Scattering_Type scattering_type = Scattering_Type::FULL);
    
private: 

    // Apply within-group and out-of-group fission
    virtual void apply_full(vector<double> &x) const override;

    // Apply only within-group fission
    virtual void apply_coherent(vector<double> &x) const override;

    void calculate_cross_sections(vector<double> &chi,
                                  vector<double> &nu,
                                  vector<double> &sigma_f) const;
};

#endif
