#ifndef Cartesian_Distance_hh
#define Cartesian_Distance_hh

#include "Distance.hh"

class Cartesian_Distance : public Distance
{
public:
    
    Cartesian_Distance(int dimension);

    virtual int dimension() const override
    {
        return dimension_;
    }
    
    virtual double distance(std::vector<double> const &r,
                            std::vector<double> const &r0) const override;
    
    virtual double d_distance(int dim,
                              std::vector<double> const &r,
                              std::vector<double> const &r0) const override;
    
    virtual double dd_distance(int dim,
                               std::vector<double> const &r,
                               std::vector<double> const &r0) const override;

    virtual std::vector<double> gradient_distance(std::vector<double> const &r,
                                                  std::vector<double> const &r0) const override;

    virtual std::vector<double> double_gradient_distance(std::vector<double> const &r,
                                                         std::vector<double> const &r0) const override;
    
    virtual double laplacian_distance(std::vector<double> const &r,
                                      std::vector<double> const &r0) const override;
    
    virtual std::string description() const override
    {
        return "cartesian";
    }

    virtual void check_class_invariants() const override;
    
private:

    int dimension_;
};

#endif
