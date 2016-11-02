#ifndef Material_hh
#define Material_hh

#include <memory>
#include <vector>

#include "pugixml.hh"

class Angular_Discretization;
class Energy_Discretization;

using std::shared_ptr;
using std::vector;

/*
  Holds cross section data for the problem
*/
class Material
{
public:
    
    Material(int index,
             shared_ptr<Angular_Discretization> angular_discretization,
             shared_ptr<Energy_Discretization> energy_discretization,
             vector<double> const &sigma_t,
             vector<double> const &sigma_s,
             vector<double> const &nu,
             vector<double> const &sigma_f,
             vector<double> const &chi,
             vector<double> const &internal_source);

    // Material index
    virtual int index() const
    {
        return index_;
    }

    // Total cross section
    virtual vector<double> const &sigma_t() const
    {
        return sigma_t_;
    }
    
    // Scattering cross section
    virtual vector<double> const &sigma_s() const
    {
        return sigma_s_;
    }

    // Average number of fission neutrons emitted
    virtual vector<double> const &nu() const
    {
        return nu_;
    }

    // Fission cross section
    virtual vector<double> const &sigma_f() const
    {
        return sigma_f_;
    }

    // Fission distribution function
    virtual vector<double> const &chi() const
    {
        return chi_;
    }

    virtual vector<double> const &internal_source() const
    {
        return internal_source_;
    }

    virtual shared_ptr<Angular_Discretization> angular_discretization() const
    {
        return angular_discretization_;
    }

    virtual shared_ptr<Energy_Discretization> energy_discretization() const
    {
        return energy_discretization_;
    }
    
    // Check class invariants
    virtual void check_class_invariants() const;

    // Output data to XML file
    virtual void output(pugi::xml_node &output_node) const;

private:
    
    int index_;
    
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
    
    vector<double> sigma_t_;
    vector<double> sigma_s_;
    vector<double> nu_;
    vector<double> sigma_f_;
    vector<double> chi_;
    vector<double> internal_source_;
};

#endif
