#ifndef Legendre_Function_hh
#define Legendre_Function_hh

#include <memory>
#include <vector>

#include "Meshless_Function.hh"

class Legendre;
class XML_Node;

class Legendre_Function : public Meshless_Function
{
public:

    Legendre_Function(int index,
                      int dimension,
                      std::vector<int> order,
                      std::vector<std::vector<double> > limits);

    // Functions shared for all dimensions
    virtual bool depends_on_neighbors() const override
    {
        return false;
    }
    virtual int index() const override
    {
        return index_;
    }
    virtual int dimension() const override
    {
        return dimension_;
    }
    virtual double radius() const override;
    virtual std::vector<double> position() const override
    {
        return position_;
    }
    virtual double shape() const override
    {
        return shape_;
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

    void initialize();
    
    int index_;
    int dimension_;
    double shape_;
    std::vector<int> order_;
    std::vector<std::vector<double> > limits_;
    std::vector<double> position_;
    std::vector<std::shared_ptr<Legendre> > functions_;
};

#endif
