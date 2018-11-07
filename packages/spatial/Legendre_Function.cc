#include "Legendre_Function.hh"

#include <limits>

#include "Distance.hh"
#include "Legendre.hh"
#include "XML_Node.hh"

using namespace std;

Legendre_Function::
Legendre_Function(int index,
                  int dimension,
                  std::vector<int> order,
                  std::vector<std::vector<double> > limits):
    index_(index),
    dimension_(dimension),
    order_(order),
    limits_(limits)
{
    initialize();
    check_class_invariants();
}

void Legendre_Function::
initialize()
{
    // Set position
    position_.resize(dimension_);
    for (int d = 0; d < dimension_; ++d)
    {
        position_[d] = 0.5 * (limits_[d][1] + limits_[d][0]);
    }
    
    // Get 1D functions
    functions_.resize(dimension_);
    for (int d = 0; d < dimension_; ++d)
    {
        functions_[d] = make_shared<Legendre>(order_[d],
                                              limits_[d]);
    }

    // Calculate pseudo-shape parameter for tau calculations by averaging side lengths
    double delta = 0;
    for (int d = 0; d < dimension_; ++d)
    {
        delta += limits_[d][1] - limits_[d][0];
    }
    shape_ = dimension_ / delta;
}

double Legendre_Function::
radius() const
{
    return 0.5 * numeric_limits<double>::max();    
}

double Legendre_Function::
value(vector<double> const &r) const
{
    double val = 1;
    for (int d = 0; d < dimension_; ++d)
    {
        val *= functions_[d]->value(r[d]);
    }

    return val;
}

double Legendre_Function::
d_value(int dim,
        vector<double> const &r) const
{
    double val = 1;
    for (int d = 0; d < dimension_; ++d)
    {
        if (d == dim)
        {
            val *= functions_[d]->d_value(r[d]);
        }
        else
        {
            val *= functions_[d]->value(r[d]);
        }
    }

    return val;
}


double Legendre_Function::
dd_value(int dim,
         vector<double> const &r) const
{
    double val = 1;
    for (int d = 0; d < dimension_; ++d)
    {
        if (d == dim)
        {
            val *= functions_[d]->dd_value(r[d]);
        }
        else
        {
            val *= functions_[d]->value(r[d]);
        }
    }

    return val;
}

vector<double> Legendre_Function::
gradient_value(vector<double> const &r) const               
{
    vector<double> res(dimension_);
    
    for (int d = 0; d < dimension_; ++d)
    {
        res[d] = d_value(d,
                         r);
    }

    return res;
}

double Legendre_Function::
laplacian_value(vector<double> const &r) const
{
    double res = 0;
    for (int d = 0; d < dimension_; ++d)
    {
        res += dd_value(d,
                        r);
    }
    
    return res;
}

void Legendre_Function::
output(XML_Node output_node) const
{
    XML_Node rbf_node = output_node.append_child("legendre_function");

    rbf_node.set_child_value(index_, "index");
    rbf_node.set_child_vector(order_, "order");
}

void Legendre_Function::
check_class_invariants() const
{
    Assert(order_.size() == dimension_);
    Assert(limits_.size() == dimension_);
    Assert(position_.size() == dimension_);
    Assert(functions_.size() == dimension_);
}

