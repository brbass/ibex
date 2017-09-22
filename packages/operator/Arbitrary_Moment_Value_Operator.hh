#ifndef Arbitrary_Moment_Value_Operator_hh
#define Arbitrary_Moment_Value_Operator_hh

#include "Vector_Operator.hh"

#include <memory>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;
class Weak_Spatial_Discretization;

class Arbitrary_Moment_Value_Operator : public Vector_Operator
{
public:

    Arbitrary_Moment_Value_Operator(std::shared_ptr<Weak_Spatial_Discretization> spatial,
                                    std::shared_ptr<Angular_Discretization> angular,
                                    std::shared_ptr<Energy_Discretization> energy,
                                    std::vector<std::vector<double> > const &evaluation_points);

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
        return "Arbitrary_Moment_Value_Operator";
    }
    
private:

    int row_size_;
    int column_size_;
    int number_of_evaluation_points_;
    std::vector<std::vector<double> > evaluation_points_;
    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
};

#endif
