#ifndef Cross_Section_hh
#define Cross_Section_hh

#include <memory>
#include <string>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;
class XML_Node;

class Cross_Section
{
public:
    
    enum class Angular_Dependence
    {
        NONE,
        SCATTERING_MOMENTS,
        MOMENTS,
        ORDINATES
    };
    
    enum class Energy_Dependence
    {
        NONE,
        GROUP,
        GROUP_TO_GROUP
    };
    
    Cross_Section(Angular_Dependence angular_dependence,
                  Energy_Dependence energy_dependence,
                  std::shared_ptr<Angular_Discretization> angular_discretization,
                  std::shared_ptr<Energy_Discretization> energy_discretization,
                  std::vector<double> const &data);

    virtual Angular_Dependence angular_dependence() const
    {
        return angular_dependence_;
    }
    virtual Energy_Dependence energy_dependence() const
    {
        return energy_dependence_;
    }
    virtual std::vector<double> const &data() const
    {
        return data_;
    }
    
    virtual int size() const;
    virtual int angular_size() const;
    virtual int energy_size() const;
    virtual std::string angular_string() const;
    virtual std::string energy_string() const;
    
    virtual void check_class_invariants() const;
    virtual void output(XML_Node output_node) const;
    
private:
    
    Angular_Dependence angular_dependence_;
    Energy_Dependence energy_dependence_;
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;
    std::vector<double> data_;
};

#endif
