#ifndef Integral_Error_Operator_hh
#define Integral_Error_Operator_hh

#include "Integration_Mesh.hh"
#include "Vector_Operator.hh"

#include <functional>
#include <memory>
#include <string>
#include <vector>

class Angular_Discretization;
template<class T1, class T2> class Conversion;
class Energy_Discretization;
class Weak_Spatial_Discretization;

class Integral_Error_Operator : public Vector_Operator
{
public:

    struct Options
    {
        enum class Norm
        {
            INTEGRAL,
            L1,
            L2,
            LINF
        };

        enum class Angular
        {
            NONE,
            MOMENTS,
            ORDINATES,
        };

        enum class Energy
        {
            NONE,
            GROUP
        };
        
        std::shared_ptr<Conversion<Norm, std::string> > norm_conversion() const;
        std::shared_ptr<Conversion<Angular, std::string> > angular_conversion() const;
        std::shared_ptr<Conversion<Energy, std::string> > energy_conversion() const;
        
        Norm norm;
        Angular angular;
        Energy energy;
    };
    
    Integral_Error_Operator(Options options,
                            std::shared_ptr<Integration_Mesh_Options> integration_options,
                            std::shared_ptr<Weak_Spatial_Discretization> spatial,
                            std::shared_ptr<Angular_Discretization> angular,
                            std::shared_ptr<Energy_Discretization> energy,
                            std::function<std::vector<double>(std::vector<double> const&)> const& solution);
    
    virtual int row_size() const override
    {
        return row_size_;
    }
    virtual int column_size() const override
    {
        return column_size_;
    }

    virtual void check_class_invariants() const override;
    virtual std::string description() const override
    {
        return "Integral_Error_Operator";
    }

    virtual int angular_size() const
    {
        return angular_size_;
    }
    
private:

    virtual void apply(std::vector<double> &x) const override;
    
    void get_flux(std::shared_ptr<Integration_Cell> const cell,
                  std::vector<double> const &b_val,
                  std::vector<double> const &coeff,
                  std::vector<double> &flux) const;
    
    int row_size_;
    int column_size_;
    int angular_size_;
    int energy_size_;
    int number_per_point_;

    Options options_;
    std::shared_ptr<Integration_Mesh_Options> integration_options_;
    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
    std::function<std::vector<double>(std::vector<double> const&)> const solution_;
    
    mutable bool initialized_;
    mutable std::shared_ptr<Integration_Mesh> mesh_;
};

#endif
   
