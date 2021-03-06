#ifndef Cross_Section_hh
#define Cross_Section_hh

#include <memory>
#include <string>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;
template<class T1, class T2> class Conversion;
class XML_Node;

class Cross_Section
{
public:

    struct Dependencies
    {
        enum class Angular
        {
            NONE,
            SCATTERING_MOMENTS,
            MOMENTS,
            ORDINATES
        };
        
        enum class Energy
        {
             NONE,
             GROUP,
             GROUP_TO_GROUP
        };
        
        enum class Dimensional
        {
            NONE,
            SUPG
        };

        enum class Spatial
        {
            WEIGHT,
            BASIS,
            BASIS_WEIGHT
        };
        
        Angular angular = Angular::NONE;
        Energy energy = Energy::NONE;
        Dimensional dimensional = Dimensional::NONE;
        Spatial spatial = Spatial::WEIGHT;

        int number_of_basis_functions = -1;
        
        std::shared_ptr<Conversion<Angular, std::string> > angular_conversion() const;
        std::shared_ptr<Conversion<Energy, std::string> > energy_conversion() const;
        std::shared_ptr<Conversion<Dimensional, std::string> > dimensional_conversion() const;
        std::shared_ptr<Conversion<Spatial, std::string> > spatial_conversion() const;
    };
    
    Cross_Section(Dependencies dependencies,
                  std::shared_ptr<Angular_Discretization> angular_discretization,
                  std::shared_ptr<Energy_Discretization> energy_discretization,
                  std::vector<double> const &data);

    virtual Dependencies dependencies() const
    {
        return dependencies_;
    }
    virtual std::vector<double> const &data() const
    {
        return data_;
    }
    
    virtual int size() const;
    virtual int angular_size() const;
    virtual int energy_size() const;
    virtual int dimensional_size() const;
    virtual int spatial_size() const;
    virtual std::string angular_string() const;
    virtual std::string energy_string() const;
    virtual std::string dimensional_string() const;
    virtual std::string spatial_string() const;
    
    virtual void check_class_invariants() const;
    virtual void output(XML_Node output_node) const;
    
private:
    
    Dependencies dependencies_;
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;
    std::vector<double> data_;
};

#endif
