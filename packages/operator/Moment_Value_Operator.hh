#ifndef Moment_Value_Operator_hh
#define Moment_Value_Operator_hh

#include "Square_Vector_Operator.hh"

#include <memory>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;
class Weak_Spatial_Discretization;

class Moment_Value_Operator : public Square_Vector_Operator
{
public:

    Moment_Value_Operator(std::shared_ptr<Weak_Spatial_Discretization> spatial,
                          std::shared_ptr<Angular_Discretization> angular,
                          std::shared_ptr<Energy_Discretization> energy,
                          bool weighted);

    virtual int size() const
    {
        return size_;
    }
    
    virtual void check_class_invariants() const override;
    virtual void apply(std::vector<double> &x) const override;

    virtual std::string description() const override
    {
        return "Moment_Value_Operator";
    }
    
private:

    bool weighted_;
    int size_;
    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
};

#endif
