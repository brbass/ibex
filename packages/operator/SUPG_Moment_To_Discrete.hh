#ifndef SUPG_Moment_To_Discrete_hh
#define SUPG_Moment_To_Discrete_hh

#include <memory>

#include "Vector_Operator.hh"

class Angular_Discretization;
class Dimensional_Moments;
class Energy_Discretization;
class Weak_Spatial_Discretization;

class SUPG_Moment_To_Discrete : public Vector_Operator
{
public:

// Constructor
    SUPG_Moment_To_Discrete(std::shared_ptr<Weak_Spatial_Discretization> spatial_discretization,
                            std::shared_ptr<Angular_Discretization> angular_discretization,
                            std::shared_ptr<Energy_Discretization> energy_discretization,
                            bool include_double_dimensional_moments);

    virtual void check_class_invariants() const override;

    virtual int row_size() const override
    {
        return row_size_;
    }
    virtual int column_size() const override
    {
        return column_size_;
    }
    virtual std::string description() const override
    {
        return "SUPG_Moment_To_Discrete";
    }
    
private:

    virtual void apply(std::vector<double> &x) const override;

    bool include_double_dimensional_moments_;
    int row_size_;
    int column_size_;
    int local_number_of_dimensional_moments_;
    std::shared_ptr<Dimensional_Moments> dimensional_moments_;
    std::shared_ptr<Weak_Spatial_Discretization> spatial_discretization_;
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;
};

#endif
