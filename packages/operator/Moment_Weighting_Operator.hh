#ifndef Moment_Weighting_Operator_hh
#define Moment_Weighting_Operator_hh

#include "Weighting_Operator.hh"

#include <memory>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;
class Weak_Spatial_Discretization;

class Moment_Weighting_Operator : public Weighting_Operator
{
public:

    Moment_Weighting_Operator(std::shared_ptr<Weak_Spatial_Discretization> spatial,
                              std::shared_ptr<Angular_Discretization> angular,
                              std::shared_ptr<Energy_Discretization> energy,
                              Options options);    
    virtual int row_size() const override
    {
        return row_size_;
    }
    virtual int column_size() const override
    {
        return column_size_;
    }
    
    virtual void check_class_invariants() const override;
    virtual void apply(std::vector<double> &x) const override;
    virtual std::string description() const override
    {
        return "Moment_Weighting_Operator";
    }
    
private:

    int row_size_;
    int column_size_;
};

#endif
