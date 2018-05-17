#ifndef Quadratic_MLS_Function_hh
#define Quadratic_MLS_Function_hh

#include <memory>
#include <vector>

#include "Meshless_Function.hh"

class Quadratic_MLS_Normalization;
class XML_Node;
template<class T> class Dense_Solver;

/*
  Moving least squares function

  neighbor_functions: the weight functions; neighbor_functions[0] should be the
  weight function for the current MLS function
*/
class Quadratic_MLS_Function : public Meshless_Function
{
public:

    Quadratic_MLS_Function(std::vector<std::shared_ptr<Meshless_Function> > neighbor_functions);

    // Meshless_Function methods
    virtual bool depends_on_neighbors() const override
    {
        return true;
    }
    virtual int index() const override
    {
        return index_;
    }
    virtual int dimension() const override
    {
        return dimension_;
    }
    virtual double radius() const override
    {
        return radius_;
    }
    virtual double shape() const override
    {
        return function_->shape();
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
    virtual void values(std::vector<double> const &r,
                        std::vector<int> &indices,
                        std::vector<double> &values) const override;
    virtual void gradient_values(std::vector<double> const &r,
                                 std::vector<int> &indices,
                                 std::vector<double> &vals,
                                 std::vector<std::vector<double> > &grad_vals) const override;
    virtual std::shared_ptr<Meshless_Function> base_function() override
    {
        return function_;
    }
    virtual std::shared_ptr<Meshless_Normalization> normalization() const override;
    virtual void output(XML_Node output_node) const override;
    virtual void check_class_invariants() const override;

private:

    // Quadratic_MLS_Function methods
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
    int index_;
    int dimension_;
    int number_of_polynomials_;
    int number_of_functions_;
    double radius_;
    std::vector<double> position_;
    std::shared_ptr<Meshless_Function> function_;
    std::vector<std::shared_ptr<Meshless_Function> > neighbor_functions_;
    std::shared_ptr<Quadratic_MLS_Normalization> normalization_;
    std::shared_ptr<Dense_Solver<double> > solver_;
};

#endif
