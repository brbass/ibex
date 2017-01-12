#include "Linear_MLS_Function.hh"

#include "Linear_Algebra.hh"
#include "Matrix_Functions.hh"
#include "Vector_Functions.hh"

namespace la = Linear_Algebra;
namespace mf = Matrix_Functions;
namespace vf = Vector_Functions;

using std::shared_ptr;
using std::vector;

Linear_MLS_Function::
Linear_MLS_Function(vector<shared_ptr<Meshless_Function> > neighbor_functions):
    dimension_(neighbor_functions[0]->dimension()),
    number_of_functions_(neighbor_functions.size()),
    position_(neighbor_functions[0]->position()),
    function_(neighbor_functions[0]),
    neighbor_functions_(neighbor_functions)
{
}

double Linear_MLS_Function::
radius() const
{
    return function_->radius();
}

double Linear_MLS_Function::
basis(vector<double> const &r) const
{
    int num_poly = dimension_ + 1;

    // Get A
    vector<double> a_mat(num_poly * num_poly);
    get_a(r,
          a_mat);

    // Get B
    vector<double> b_vec(num_poly);
    get_b(r,
          b_vec);
    
    // Get x = A^-1 B
    vector<double> x_vec(num_poly);
    la::linear_solve(a_mat,
                     b_vec,
                     x_vec);
    
    // Get p
    vector<double> p_vec(num_poly);
    get_polynomial(r,
                   p_vec);

    // Return result
    return vf::dot(p_vec, x_vec);
}

double Linear_MLS_Function::
d_basis(int dim,
        vector<double> const &r) const
{
    // Get A
    vector<double> a_mat(num_poly * num_poly);
    vector<double> d_a_mat(num_poly * num_poly);
    get_d_a(dim,
            r,
            a_mat,
            d_a_mat);
    
    // Get A inverse
    vector<double> a_inv_mat(num_poly * num_poly);
    la::direct_inverse(a_mat,
                       a_inv_mat);
    vector<double> d_a_inv_mat(num_poly * num_poly);
    d_a_inv_mat = mf::square_matrix_matrix_product(num_poly,
                                                   d_a_mat,
                                                   a_inv_mat);
    d_a_inv_mat = mf::square_matrix_matrix_product(num_poly,
                                                   a_inv_mat,
                                                   d_a_inv_mat);
    d_a_inv_mat = vf::multiply(d_a_inv_mat,
                               -1.);
    
    // Get B
    vector<double> b_vec(num_poly);
    vector<double> d_b_vec(num_poly);
    get_d_b(dim,
            r,
            b_vec,
            d_b_vec);
    
    // Get p
    vector<double> p_vec(num_poly);
    get_polynomial(r,
                   p_vec);
    vector<double> d_p_vec(num_poly);
    get_d_polynomial(r,
                     p_vec);

    // Calculate terms
    double t1 = vf::dot(d_p_vec,
                        mf::square_matrix_vector_product(num_poly,
                                                         a_inv_mat,
                                                         b_vec));
    double t2 = vf::dot(p_vec,
                        mf::square_matrix_vector_product(num_poly,
                                                         d_a_inv_mat,
                                                         b_vec));
    double t3 = vf::dot(p_vec,
                        mf::square_matrix_vector_product(num_poly,
                                                         a_inv_mat,
                                                         d_b_vec));

    return t1 + t2 + t3;
}


double Linear_MLS_Function::
dd_basis(int dim,
         vector<double> const &r) const
{
    AssertMsg(false, "dd_basis not available for Linear_MLS");
    return -1.;
}

vector<double> Linear_MLS_Function::
gradient_basis(vector<double> const &r) const               
{
    vector<double> result(dimension_);

    for (int d = 0; d < dimension_; ++d)
    {
        result[d] = d_basis(d,
                            r);
    }

    return result;
}

double Linear_MLS_Function::
laplacian(vector<double> const &r)
{
    AssertMsg(false, "laplacian not available for Linear_MLS");
}

void Linear_MLS_Function::
output(XML_Node output_node) const
{
    
}

void Linear_MLS_Function::
check_class_invariants() const
{
    Assert(dimension_ >= 1);
    Assert(number_of_functions_ >= 3);
    Assert(function_);
    for (shared_ptr<Meshless_Function> func : neighbor_functions_)
    {
        Assert(func);
    }
}

void Linear_MLS_Function::
get_polynomial(vector<double> const &position,
               vector<double> &poly) const
{
    poly.resize(dimension_ + 1);
    
    poly[0] = 1.;
    for (int d = 0; d < dimension_; ++d)
    {
        poly[d+1] = position[d];
    }
}

void Linear_MLS_Function::
get_d_polynomial(int dim,
                 vector<double> const &position,
                 vector<double> &d_poly) const
{
    d_poly.assign(dimension_ + 1, 0);
    
    d_poly[dim] = 1.;
}

void Linear_MLS_Function::
get_a(vector<double> const &position,
      vector<double> &a) const
{
    int num_poly = dimension_ + 1;

    a.assign(num_poly * num_poly, 0);
    for (int i = 0; i < number_of_functions_; ++i)
    {
        shared_ptr<Meshless_Function> func = neighbor_functions_[i];
        vector<double> poly;
        get_polynomial(func->position(),
                       poly);
        double weight = func->basis(position);
        
        for (int j = 0; j < num_poly; ++j)
        {
            for (int k = 0; k < num_poly; ++k)
            {
                int l = j + num_poly * k;
                double polyprod = poly[j] * poly[k];
                a[l] += weight * polyprod;
            }
        }
    }
}

void Linear_MLS_Function::
get_d_a(int dim,
        vector<double> const &position,
        vector<double> &a,
        vector<double> &d_a) const
{
    int num_poly = dimension_ + 1;
    
    a.assign(num_poly * num_poly, 0);
    d_a.assign(num_poly * num_poly, 0);
    
    for (int i = 0; i < number_of_functions_; ++i)
    {
        shared_ptr<Meshless_Function> func = neighbor_functions_[i];
        vector<double> poly;
        get_polynomial(func->position(),
                       poly);
        double weight = func->basis(position);
        double d_weight = func->d_basis(dim,
                                        position);
        
        for (int j = 0; j < num_poly; ++j)
        {
            for (int k = 0; k < num_poly; ++k)
            {
                int l = j + num_poly * k;
                double polyprod = poly[j] * poly[k];
                a[l] += weight * polyprod;
                d_a[l] += d_weight * polyprod;
            }
        }
    }
}

void Linear_MLS_Function::
get_b(vector<double> const &position,
      vector<double> &b) const
{
    b.resize(num_poly);
    
    vector<double> poly;
    get_polynomial(function_->position,
                   poly);
    double weight = function_->basis(position);
    
    for (int i = 0; i < num_poly; ++i)
    {
        b[i] = poly[i] * weight;
    }
}

void Linear_MLS_Function::
get_d_b(int dim,
        vector<double> const &position,
        vector<double> &b,
        vector<double> &d_b) const
{
    b.resize(num_poly);
    d_b.resize(num_poly);
    
    vector<double> poly;
    get_polynomial(function_->position,
                   poly);
    double weight = function_->basis(position);
    double d_weight = function_->d_basis(dim,
                                         position);
    
    for (int i = 0; i < num_poly; ++i)
    {
        b[i] = poly[i] * weight;
        d_b[i] = poly[i] * d_weight;
    }
}