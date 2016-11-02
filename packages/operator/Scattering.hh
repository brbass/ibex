#ifndef Scattering_hh
#define Scattering_hh

#include <memory>
#include <vector>

#include "Scattering_Operator.hh"

using std::shared_ptr;
using std::vector;

/*
  Applies scattering to a moment representation of the flux
*/
class Scattering : public Scattering_Operator
{
public:

    // Constructor
    Scattering(shared_ptr<Spatial_Discretization> spatial_discretization,
               shared_ptr<Angular_Discretization> angular_discretization,
               shared_ptr<Energy_Discretization> energy_discretization,
               Scattering_Type scattering_type = Scattering_Type::FULL);
    
private: 

    // Apply within-group and out-of-group scattering
    virtual void apply_full(vector<double> &x) const override;

    // Apply only within-group scattering
    virtual void apply_coherent(vector<double> &x) const override;
};

#endif
