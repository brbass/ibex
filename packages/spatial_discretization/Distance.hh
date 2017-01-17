#ifndef Distance_hh
#define Distance_hh

#include <string>
#include <vector>

class Distance
{
public:
    
    Distance();

    // Methods to be implemented by child classes
    
    virtual int dimension() const = 0;
    
    virtual double distance(std::vector<double> const &r,
                            std::vector<double> const &r0) const = 0;
    
    virtual double d_distance(int dim,
                              std::vector<double> const &r,
                              std::vector<double> const &r0) const = 0;
    
    virtual double dd_distance(int dim,
                               std::vector<double> const &r,
                               std::vector<double> const &r0) const = 0;

    virtual std::vector<double> gradient_distance(std::vector<double> const &r,
                                             std::vector<double> const &r0) const = 0;

    virtual std::vector<double> double_gradient_distance(std::vector<double> const &r,
                                                    std::vector<double> const &r0) const = 0;
    
    virtual double laplacian_distance(std::vector<double> const &r,
                                      std::vector<double> const &r0) const = 0;

    virtual std::string description() const = 0;
    
    virtual void check_class_invariants() const = 0;

    virtual bool derivative_available(int derivative) const
    {
        return derivative < 3 ? true : false;
    }
};

#endif
