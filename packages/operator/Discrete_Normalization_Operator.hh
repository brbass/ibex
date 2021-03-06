#ifndef Discrete_Normalization_Operator_hh
#define Discrete_Normalization_Operator_hh

#include "Weighting_Operator.hh"

#include <memory>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;
class Weak_Spatial_Discretization;

class Discrete_Normalization_Operator : public Weighting_Operator
{
public:

    Discrete_Normalization_Operator(std::shared_ptr<Weak_Spatial_Discretization> spatial,
                                    std::shared_ptr<Angular_Discretization> angular,
                                    std::shared_ptr<Energy_Discretization> energy,
                                    Options options = Options());
    
    virtual int row_size() const override
    {
        return size_;
    }
    virtual int column_size() const override
    {
        return size_;
    }
    
    virtual void check_class_invariants() const override;
    virtual void apply(std::vector<double> &x) const override;
    virtual std::string description() const override
    {
        return "Discrete_Normalization_Operator";
    }
    
private:

    int size_;
};

#endif
