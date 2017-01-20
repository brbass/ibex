#ifndef Material_hh
#define Material_hh

#include <memory>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;
class Cross_Section;
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
             std::shared_ptr<Cross_Section> sigma_t,
             std::shared_ptr<Cross_Section> sigma_s,
             std::shared_ptr<Cross_Section> nu,
             std::shared_ptr<Cross_Section> sigma_f,
             std::shared_ptr<Cross_Section> chi,
             std::shared_ptr<Cross_Section> internal_source);
    
    // Material index
    virtual int index() const
    {
        return index_;
    }

    // Total cross section
    virtual std::shared_ptr<Cross_Section> sigma_t() const
    {
        return sigma_t_;
    }
    
    // Scattering cross section
    virtual std::shared_ptr<Cross_Section> sigma_s() const
    {
        return sigma_s_;
    }

    // Average number of fission neutrons emitted
    virtual std::shared_ptr<Cross_Section> nu() const
    {
        return nu_;
    }

    // Fission cross section
    virtual std::shared_ptr<Cross_Section> sigma_f() const
    {
        return sigma_f_;
    }

    // Fission distribution function
    virtual std::shared_ptr<Cross_Section> chi() const
    {
        return chi_;
    }

    // Internal source
    virtual std::shared_ptr<Cross_Section> internal_source() const
    {
        return internal_source_;
    }

    // Angular discretization
    virtual std::shared_ptr<Angular_Discretization> angular_discretization() const
    {
        return angular_discretization_;
    }

    // Energy discretization
    virtual std::shared_ptr<Energy_Discretization> energy_discretization() const
    {
        return energy_discretization_;
    }
    
    // Check class invariants
    virtual void check_class_invariants() const;

    // Output data to XML file
    virtual void output(XML_Node output_node) const;

protected:
    
    int index_;
    
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;
    
    std::shared_ptr<Cross_Section> sigma_t_;
    std::shared_ptr<Cross_Section> sigma_s_;
    std::shared_ptr<Cross_Section> nu_;
    std::shared_ptr<Cross_Section> sigma_f_;
    std::shared_ptr<Cross_Section> chi_;
    std::shared_ptr<Cross_Section> internal_source_;
};

#endif
