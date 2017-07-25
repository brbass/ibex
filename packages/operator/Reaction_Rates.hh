#ifndef Reaction_Rates_hh
#define Reaction_Rates_hh

#include <memory>

#include "Vector_Operator.hh"

class Angular_Discretization;
class Energy_Discretization;
class Weak_Spatial_Discretization;

class Reaction_Rates : public Vector_Operator
{
public:

    Reaction_Rates(std::shared_ptr<Weak_Spatial_Discretization> spatial,
                   std::shared_ptr<Angular_Discretization> angular,
                   std::shared_ptr<Energy_Discretization> energy);
    
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
        return "Reaction_Rates";
    }
    
private:
    
    virtual void apply(std::vector<double> &x) const override;

    void get_coefficients(int i,
                          std::vector<double> &sigma_t,
                          std::vector<double> &sigma_s,
                          std::vector<double> &sigma_f) const;
    
    int row_size_;
    int column_size_;
    std::shared_ptr<Weak_Spatial_Discretization> spatial_;
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
};

#endif
