#ifndef Manufactured_Integral_Operator_hh
#define Manufactured_Integral_Operator_hh

#include "Integration_Mesh.hh"
#include "Vector_Operator.hh"

#include <memory>
#include <string>
#include <vector>

class Angular_Discretization;
template<class T1, class T2> class Conversion;
class Energy_Discretization;
class Manufactured_Solution;
class Weak_Spatial_Discretization;

class Manufactured_Integral_Operator : public Vector_Operator
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

        std::shared_ptr<Conversion<Norm, std::string> > norm_conversion() const;
        
        Norm norm;
    };
    
    Manufactured_Integral_Operator(Options options,
                                   std::shared_ptr<Integration_Mesh_Options> integration_options,
                                   std::shared_ptr<Weak_Spatial_Discretization> spatial,
                                   std::shared_ptr<Angular_Discretization> angular,
                                   std::shared_ptr<Energy_Discretization> energy,
                                   std::shared_ptr<Manufactured_Solution> solution);

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
        return "Manufactured_Integral_Operator";
    }

private:

    virtual void apply(std::vector<double> &x) const override;
    
    void get_flux(std::shared_ptr<Integration_Mesh::Cell> const cell,
                  std::vector<double> const &b_val,
                  std::vector<double> const &coeff,
                  std::vector<double> &flux) const;
    
    int row_size_;
    int column_size_;

    Options options_;
    std::shared_ptr<Integration_Mesh_Options> integration_options_;
    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
    std::shared_ptr<Manufactured_Solution> solution_;
    
    mutable bool initialized_;
    mutable std::shared_ptr<Integration_Mesh> mesh_;
};

#endif
   
