#ifndef Material_hh
#define Material_hh

#include <memory>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;
class XML_Node;

/*
  Holds cross section data for the problem
*/
class Material
{
public:
    
    Material(int index,
             std::shared_ptr<Angular_Discretization> angular_discretization,
             std::shared_ptr<Energy_Discretization> energy_discretization,
             std::vector<double> const &sigma_t,
             std::vector<double> const &sigma_s,
             std::vector<double> const &nu,
             std::vector<double> const &sigma_f,
             std::vector<double> const &chi,
             std::vector<double> const &internal_source);

    // Material index
    virtual int index() const
    {
        return index_;
    }

    // Total cross section
    virtual std::vector<double> const &sigma_t() const
    {
        return sigma_t_;
    }
    
    // Scattering cross section
    virtual std::vector<double> const &sigma_s() const
    {
        return sigma_s_;
    }

    // Average number of fission neutrons emitted
    virtual std::vector<double> const &nu() const
    {
        return nu_;
    }

    // Fission cross section
    virtual std::vector<double> const &sigma_f() const
    {
        return sigma_f_;
    }

    // Fission distribution function
    virtual std::vector<double> const &chi() const
    {
        return chi_;
    }

    virtual std::vector<double> const &internal_source() const
    {
        return internal_source_;
    }

    virtual std::shared_ptr<Angular_Discretization> angular_discretization() const
    {
        return angular_discretization_;
    }

    virtual std::shared_ptr<Energy_Discretization> energy_discretization() const
    {
        return energy_discretization_;
    }
    
    // Check class invariants
    virtual void check_class_invariants() const;

    // Output data to XML file
    virtual void output(XML_Node output_node) const;

private:
    
    int index_;
    
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;
    
    std::vector<double> sigma_t_;
    std::vector<double> sigma_s_;
    std::vector<double> nu_;
    std::vector<double> sigma_f_;
    std::vector<double> chi_;
    std::vector<double> internal_source_;
};

#endif
