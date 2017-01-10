#ifndef RBF_Function_hh
#define RBF_Function_hh

#include <memory>
#include <vector>

#include "XML_Node.hh"

class Distance;
class RBF;

class RBF_Function : public Meshless_Function
{
public:

    RBF_Function(double shape,
                 std::vector<double> const &position,
                 std::shared_ptr<RBF> rbf,
                 std::shared_ptr<Distance> distance);
    
    virtual double radius() const override;
    virtual std::vector<double> position() const override;
    
    virtual double basis(std::vector<double> const &r) const override;
    virtual double d_basis(int dim,
                           std::vector<double> const &r) const override;
    virtual double dd_basis(int dim,
                            std::vector<double> const &r) const override;
    virtual std::vector<double> gradient_basis(std::vector<double> const &r) const override;
    virtual double laplacian(std::vector<double> const &r) const override;
    
    virtual void output(XML_Node output_node) const override;
    virtual void check_class_invariants() const override;
    
protected:

    double shape_;
    std::vector<double> position_;
    std::shared_ptr<RBF> rbf_;
    std::shared_ptr<Distance> distance_;
};

#endif
