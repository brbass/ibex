#ifndef Weighting_Operator_hh
#define Weighting_Operator_hh

#include "Square_Vector_Operator.hh"

#include <memory>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;
class Weak_Spatial_Discretization;

class Weighting_Operator : public Square_Vector_Operator
{
public:

    struct Options
    {
        enum class Normalized
        {
            AUTO, // Get from weight function options
            TRUE, // Do not include normalization
            FALSE // Include the normalization
        };
        
        enum class Include_SUPG
        {
            AUTO, // Get from weight function options
            TRUE, // Include SUPG terms
            FALSE // Do not include SUPG terms
        };
        
        Normalized normalized = Normalized::AUTO;
        Include_SUPG include_supg = Include_SUPG::AUTO;
    };
    
    Weighting_Operator(std::shared_ptr<Weak_Spatial_Discretization> spatial,
                       std::shared_ptr<Angular_Discretization> angular,
                       std::shared_ptr<Energy_Discretization> energy,
                       Options options);
    
    virtual void check_class_invariants() const = 0;
    virtual void apply(std::vector<double> &x) const = 0;
    
    virtual int size() const = 0;
    
protected:
    
    Options options_;
    
    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
};

#endif
