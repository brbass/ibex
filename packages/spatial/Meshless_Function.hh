#ifndef Meshless_Function_hh
#define Meshless_Function_hh

#include <vector>

class XML_Node;

class Meshless_Function
{
public:
    
    Meshless_Function();
    
    // Does meshless function depend on the values of all other functions
    virtual bool depends_on_neighbors() const = 0;

    // Index of function
    virtual int index() const = 0;
    
    // Spatial dimension
    virtual int dimension() const = 0;

    // Support of meshless function
    virtual double radius() const = 0;

    // Shape parameter
    virtual double shape() const = 0;

    // Center position
    virtual std::vector<double> position() const = 0;

    // Values and derivatives of the meshless function
    virtual double value(std::vector<double> const &r) const = 0;
    virtual double d_value(int dim,
                           std::vector<double> const &r) const = 0;
    virtual double dd_value(int dim,
                            std::vector<double> const &r) const = 0;
    virtual std::vector<double> gradient_value(std::vector<double> const &r) const = 0;
    virtual double laplacian_value(std::vector<double> const &r) const = 0;

    // Values and derivatives of all the meshless functions (including neighbors): disabled if !depends_on_neighbors
    virtual void values(std::vector<double> const &r,
                        std::vector<int> &indices,
                        std::vector<double> &values) const;
    virtual void gradient_values(std::vector<double> const &r,
                                 std::vector<int> &indices,
                                 std::vector<double> &vals,
                                 std::vector<std::vector<double> > &grad_vals) const;
    
    // Input checking and output methods
    virtual void output(XML_Node output_node) const = 0;
    virtual void check_class_invariants() const = 0;

    // Check whether point is inside support of function
    virtual bool inside_radius(std::vector<double> const &r) const;
};

#endif
