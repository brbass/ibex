#ifndef RBF_Function_hh
#define RBF_Function_hh

#include <memory>
#include <vector>

#include "Meshless_Function.hh"

class Distance;
class RBF;
class XML_Node;

class RBF_Function : public Meshless_Function
{
public:

    RBF_Function(int index,
                 double shape,
                 std::vector<double> const &position,
                 std::shared_ptr<RBF> rbf,
                 std::shared_ptr<Distance> distance);

    virtual bool depends_on_neighbors() const override
    {
        return false;
    }
    virtual int index() const override
    {
        return index_;
    }
    virtual int dimension() const override;
    virtual double radius() const override;
    virtual double shape() const override
    {
        return shape_;
    }
    virtual std::vector<double> position() const override
    {
        return position_;
    }
    virtual double value(std::vector<double> const &r) const override;
    virtual double d_value(int dim,
                           std::vector<double> const &r) const override;
    virtual double dd_value(int dim,
                            std::vector<double> const &r) const override;
    virtual std::vector<double> gradient_value(std::vector<double> const &r) const override;
    virtual double laplacian_value(std::vector<double> const &r) const override;
    
    virtual void output(XML_Node output_node) const override;
    virtual void check_class_invariants() const override;
    
protected:

    int index_;
    double shape_;
    std::vector<double> position_;
    std::shared_ptr<RBF> rbf_;
    std::shared_ptr<Distance> distance_;
};

#endif
