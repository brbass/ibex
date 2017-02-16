#ifndef Linear_MLS_Function_hh
#define Linear_MLS_Function_hh

#include <memory>
#include <vector>

#include "Meshless_Function.hh"

class XML_Node;

/*
  Moving least squares function

  neighbor_functions: the weight functions; neighbor_functions[0] should be the
  weight function for the current MLS function
*/
class Linear_MLS_Function : public Meshless_Function
{
public:

    Linear_MLS_Function(std::vector<std::shared_ptr<Meshless_Function> > neighbor_functions);

    // Meshless_Function methods
    virtual int dimension() const override
    {
        return dimension_;
    }
    virtual double radius() const override
    {
        return radius_;
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

private:

    // Linear_MLS_Function methods
    virtual void get_polynomial(std::vector<double> const &position,
                                std::vector<double> &poly) const;
    virtual void get_d_polynomial(int dim,
                                  std::vector<double> const &position,
                                  std::vector<double> &d_poly) const;
    virtual void get_a(std::vector<double> const &position,
                       std::vector<double> &a) const;
    virtual void get_d_a(int dim,
                         std::vector<double> const &position,
                         std::vector<double> &a,
                         std::vector<double> &d_a) const;
    virtual void get_b(std::vector<double> const &position,
                       std::vector<double> &b) const;
    virtual void get_d_b(int dim,
                         std::vector<double> const &position,
                         std::vector<double> &b,
                         std::vector<double> &d_b) const;
    
    // Data
    int dimension_;
    int number_of_polynomials_;
    int number_of_functions_;
    double radius_;
    std::vector<double> position_;
    std::shared_ptr<Meshless_Function> function_;
    std::vector<std::shared_ptr<Meshless_Function> > neighbor_functions_;
};

#endif
