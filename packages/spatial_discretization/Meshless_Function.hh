#ifndef Meshless_Function_hh
#define Meshless_Function_hh

#include <vector>

class XML_Node;

class Meshless_Function
{
public:

    Meshless_Function();

    virtual int dimension() const = 0;
    virtual double radius() const = 0;
    virtual std::vector<double> position() const = 0;
    virtual double basis(std::vector<double> const &r) const = 0;
    virtual double d_basis(int dim,
                           std::vector<double> const &r) const = 0;
    virtual double dd_basis(int dim,
                            std::vector<double> const &r) const = 0;
    virtual std::vector<double> gradient_basis(std::vector<double> const &r) const = 0;
    virtual double laplacian(std::vector<double> const &r) const = 0;
    virtual void output(XML_Node output_node) const = 0;
    virtual void check_class_invariants() const = 0;
}

#endif
