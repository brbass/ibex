#ifndef RBF_Function_hh
#define RBF_Function_hh

#include <memory>
#include <vector>

#include "pugixml.hh"

#include "Distance.hh"
#include "RBF.hh"

using std::shared_ptr;
using std::vector;

class RBF_Function
{
public:

    RBF_Function(shared_ptr<RBF> rbf,
                 shared_ptr<Distance> distance);
    
    virtual shared_ptr<RBF> rbf() const
    {
        return rbf_;
    }

    virtual shared_ptr<Distance> distance() const
    {
        return distance_;
    }

    virtual double basis(double shape,
                         vector<double> const &r,
                         vector<double> const &r0) const;
    
    virtual double d_basis(int dim,
                           double shape,
                           vector<double> const &r,
                           vector<double> const &r0) const;
    
    virtual double dd_basis(int dim,
                            double shape,
                            vector<double> const &r,
                            vector<double> const &r0) const;
    
    virtual vector<double> gradient_basis(double shape,
                                          vector<double> const &r,
                                          vector<double> const &r0) const;
    
    virtual double laplacian(double shape,
                             vector<double> const &r,
                             vector<double> const &r0) const;
    
    virtual void output(pugi::xml_node &output_node) const;

    virtual void check_class_invariants() const;

    virtual bool derivative_available(int derivative) const
    {
        return (derivative < 3 ? true : false
                && rbf_->derivative_available(derivative)
                && distance_->derivative_available(derivative));
    }
    
protected:

    shared_ptr<RBF> rbf_;
    shared_ptr<Distance> distance_;
};

#endif
