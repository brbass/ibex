#ifndef Weighting_Operator_hh
#define Weighting_Operator_hh

#include "Square_Vector_Operator.hh"

#include <memory>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;
template<class T1, class T2> class Conversion;
class Weak_Spatial_Discretization;

class Weighting_Operator : public Vector_Operator
{
public:
    
    struct Options
    {
        enum class Normalization
        {
            AUTO, // Get from weight function options
            TRUE, // Include normalization
            FALSE // Do not include normalization
        };
        
        enum class Include_SUPG
        {
            AUTO, // Get from weight function options
            TRUE, // Include SUPG terms
            FALSE // Do not include SUPG terms
        };
        
        Normalization normalization = Normalization::AUTO;
        Include_SUPG include_supg = Include_SUPG::AUTO;
        
        std::shared_ptr<Conversion<Normalization, std::string> > normalization_conversion() const;
        std::shared_ptr<Conversion<Include_SUPG, std::string> > include_supg_conversion() const;
    };
    
    Weighting_Operator(std::shared_ptr<Weak_Spatial_Discretization> spatial,
                       std::shared_ptr<Angular_Discretization> angular,
                       std::shared_ptr<Energy_Discretization> energy,
                       Options options);
    
    virtual void check_class_invariants() const override = 0;
    
    virtual int row_size() const override = 0;
    virtual int column_size() const override = 0;
    
    virtual std::string description() const override = 0;
    
protected:
    
    Options options_;

    int local_number_of_dimensional_moments_;
    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;

private:
    
    virtual void apply(std::vector<double> &x) const override = 0;
    
};

#endif
