#ifndef Scattering_Operator_hh
#define Scattering_Operator_hh

#include <memory>
#include <vector>

#include "Square_Vector_Operator.hh"

class Angular_Discretization;
class Energy_Discretization;
template<class T1, class T2> class Conversion;
class Spatial_Discretization;

/*
  Pure virtual class to apply scattering to a moment representation of the flux
*/
class Scattering_Operator : public Square_Vector_Operator
{
public:

    // Types of scattering
    struct Options
    {
        enum class Scattering_Type
        {
            COHERENT, // within-group scattering
            INCOHERENT, // out-of-group scattering
            FULL // within-group and out-of-group scattering
        };
        
        Scattering_Type scattering_type = Scattering_Type::FULL;

        std::shared_ptr<Conversion<Scattering_Type, std::string> > scattering_type_conversion() const;
    };

    // Constructor
    Scattering_Operator(std::shared_ptr<Spatial_Discretization> spatial_discretization,
                        std::shared_ptr<Angular_Discretization> angular_discretization,
                        std::shared_ptr<Energy_Discretization> energy_discretization,
                        Options options);
    
    virtual int size() const override
    {
        return size_;
    }
    
    virtual void check_class_invariants() const override = 0;;
    virtual std::string description() const override = 0;
    
protected:

    int size_;
    
    // Type of scattering
    Options options_;

    std::shared_ptr<Spatial_Discretization> spatial_discretization_;
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;

private: 
    
    // Apply scattering of chosen type
    virtual void apply(std::vector<double> &x) const override;

    // Apply within-group and out-of-group scattering
    virtual void apply_full(std::vector<double> &x) const = 0;
    
    // Apply only within-group scattering
    virtual void apply_coherent(std::vector<double> &x) const = 0;
    
    // Apply only out-of-group scattering
    virtual void apply_incoherent(std::vector<double> &x) const;
};

#endif
