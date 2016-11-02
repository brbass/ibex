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
    
    virtual double distance(vector<double> const &r,
                            vector<double> const &r0) const override;
    
    virtual double d_distance(int dim,
                              vector<double> const &r,
                              vector<double> const &r0) const override;
    
    virtual double dd_distance(int dim,
                               vector<double> const &r,
                               vector<double> const &r0) const override;

    virtual vector<double> gradient_distance(vector<double> const &r,
                                             vector<double> const &r0) const override;

    virtual vector<double> double_gradient_distance(vector<double> const &r,
                                                    vector<double> const &r0) const override;
    
    virtual double laplacian_distance(vector<double> const &r,
                                      vector<double> const &r0) const override;
    
    virtual string description() const override
    {
        return "cartesian";
    }

    virtual void check_class_invariants() const override;
    
private:

    int dimension_;
};

#endif
