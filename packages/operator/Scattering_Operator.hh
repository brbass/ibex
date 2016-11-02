#ifndef Scattering_Operator_hh
#define Scattering_Operator_hh

#include <memory>
#include <vector>

#include "Vector_Operator.hh"

class Angular_Discretization;
class Energy_Discretization;
class Spatial_Discretization;

using std::shared_ptr;
using std::vector;

/*
  Pure virtual class to apply scattering to a moment representation of the flux
*/
class Scattering_Operator : public Vector_Operator
{
public:

    // Types of scattering
    enum class Scattering_Type
    {
        COHERENT,
        INCOHERENT,
        FULL
    };

    // Constructor
    Scattering_Operator(shared_ptr<Spatial_Discretization> spatial_discretization,
                        shared_ptr<Angular_Discretization> angular_discretization,
                        shared_ptr<Energy_Discretization> energy_discretization,
                        Scattering_Type scattering_type = Scattering_Type::FULL);

    virtual void check_class_invariants() const override;
    
protected:

    // Type of scattering
    Scattering_Type scattering_type_;

    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;

private: 

    // Apply scattering of chosen type
    virtual void apply(vector<double> &x) const override;

    // Apply within-group and out-of-group scattering
    virtual void apply_full(vector<double> &x) const = 0;
    
    // Apply only within-group scattering
    virtual void apply_coherent(vector<double> &x) const = 0;
    
    // Apply only out-of-group scattering
    virtual void apply_incoherent(vector<double> &x) const;
};

#endif
