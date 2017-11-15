#include "Quadratic_MLS_Function.hh"

#include "Dense_Solver_Factory.hh"
#include "Quadratic_MLS_Normalization.hh"
#include "Matrix_Functions.hh"
#include "Vector_Functions.hh"
#include "XML_Node.hh"

namespace mf = Matrix_Functions;
namespace vf = Vector_Functions;

using std::make_shared;
using std::shared_ptr;
using std::vector;

Quadratic_MLS_Function::
Quadratic_MLS_Function(vector<shared_ptr<Meshless_Function> > neighbor_functions):
    index_(neighbor_functions[0]->index()),
    dimension_(neighbor_functions[0]->dimension()),
    number_of_functions_(neighbor_functions.size()),
    position_(neighbor_functions[0]->position()),
    function_(neighbor_functions[0]),
    neighbor_functions_(neighbor_functions)
{
    switch (dimension_)
    {
    case 1:
        number_of_polynomials_ = 3;
        break;
    case 2:
        number_of_polynomials_ = 6;
        break;
    case 3:
        number_of_polynomials_ = 10;
        break;
    default:
        AssertMsg(false, "dimension incorrect");
    }
    radius_ = function_->radius();
    normalization_ = make_shared<Quadratic_MLS_Normalization>(dimension_);

    Dense_Solver_Factory factory;
    solver_ = factory.get_solver(number_of_polynomials_);
    
    check_class_invariants();
}

double Quadratic_MLS_Function::
value(vector<double> const &r) const
{
    if (!function_->inside_radius(r))
    {
        return 0.;
    }

    // Get A
    vector<double> a_mat(number_of_polynomials_ * number_of_polynomials_);
    get_a(r,
          a_mat);
    
    // Get B
    vector<double> b_vec(number_of_polynomials_);
    get_b(r,
          b_vec);
    
    // Get x = A^-1 B
    vector<double> x_vec(number_of_polynomials_);
    solver_->solve(a_mat,
                   b_vec,
                   x_vec);
    
    // Get p
    vector<double> p_vec(number_of_polynomials_);
    get_polynomial(r,
                   p_vec);

    // Return result
    return vf::dot(p_vec, x_vec);
}

double Quadratic_MLS_Function::
d_value(int dim,
        vector<double> const &r) const
{
    // Check if value is inside basis radius
    if (!function_->inside_radius(r))
    {
        return 0.;
    }
    
    // Get A
    vector<double> a_mat(number_of_polynomials_ * number_of_polynomials_);
    vector<double> d_a_mat(number_of_polynomials_ * number_of_polynomials_);
    get_d_a(dim,
            r,
            a_mat,
            d_a_mat);

    // Get A inverse
    vector<double> a_inv_mat(number_of_polynomials_ * number_of_polynomials_);
    solver_->inverse(a_mat,
                     a_inv_mat);
    vector<double> d_a_inv_mat(number_of_polynomials_ * number_of_polynomials_);
    d_a_inv_mat = mf::square_matrix_matrix_product(number_of_polynomials_,
                                                   d_a_mat,
                                                   a_inv_mat);
    d_a_inv_mat = mf::square_matrix_matrix_product(number_of_polynomials_,
                                                   a_inv_mat,
                                                   d_a_inv_mat);
    d_a_inv_mat = vf::multiply(d_a_inv_mat,
                               -1.);
    
    // Get B
    vector<double> b_vec(number_of_polynomials_);
    vector<double> d_b_vec(number_of_polynomials_);
    get_d_b(dim,
            r,
            b_vec,
            d_b_vec);
    
    // Get p
    vector<double> p_vec(number_of_polynomials_);
    get_polynomial(r,
                   p_vec);
    vector<double> d_p_vec(number_of_polynomials_);
    get_d_polynomial(dim,
                     r,
                     d_p_vec);
    
    // Calculate terms
    double t1 = vf::dot(d_p_vec,
                        mf::square_matrix_vector_product(number_of_polynomials_,
                                                         a_inv_mat,
                                                         b_vec));
    double t2 = vf::dot(p_vec,
                        mf::square_matrix_vector_product(number_of_polynomials_,
                                                         d_a_inv_mat,
                                                         b_vec));
    double t3 = vf::dot(p_vec,
                        mf::square_matrix_vector_product(number_of_polynomials_,
                                                         a_inv_mat,
                                                         d_b_vec));
    return t1 + t2 + t3;
}


double Quadratic_MLS_Function::
dd_value(int dim,
         vector<double> const &r) const
{
    AssertMsg(false, "dd_value not available for Quadratic_MLS");
    return -1.;
}

vector<double> Quadratic_MLS_Function::
gradient_value(vector<double> const &r) const               
{
    // Check if value is inside basis radius
    if (!function_->inside_radius(r))
    {
        return vector<double>(dimension_, 0.);
    }
    
    // Get dimensional derivatives separately
    vector<double> result(dimension_);

    for (int d = 0; d < dimension_; ++d)
    {
        result[d] = d_value(d,
                            r);
    }

    return result;
}

double Quadratic_MLS_Function::
laplacian_value(vector<double> const &r) const
{
    AssertMsg(false, "laplacian not available for Quadratic_MLS");
    return -1.;
}

void Quadratic_MLS_Function::
values(vector<double> const &position,
       vector<int> &indices,
       vector<double> &vals) const
{
    // Check if position is inside weight function radius
    Assert(function_->inside_radius(position));
    
    // Resize indices and values
    indices.resize(number_of_functions_);
    vals.resize(number_of_functions_);
    
    // Get indices of basis functions
    for (int i = 0; i < number_of_functions_; ++i)
    {
        indices[i] = neighbor_functions_[i]->index();
    }
    
    // Get meshless function values
    vector<double> base_values(number_of_functions_);
    for (int i = 0; i < number_of_functions_; ++i)
    {
        base_values[i] = neighbor_functions_[i]->value(position);
    }

    // Get center positions
    vector<vector<double> > center_positions(number_of_functions_);
    for (int i = 0; i < number_of_functions_; ++i)
    {
        center_positions[i] = neighbor_functions_[i]->position();
    }

    // Get values
    normalization_->get_values(position,
                               center_positions,
                               base_values,
                               vals);
}

void Quadratic_MLS_Function::
gradient_values(vector<double> const &position,
                vector<int> &indices,
                vector<double> &vals,
                vector<vector<double> > &grad_vals) const
{
    // Check if position is inside weight function radius
    Assert(function_->inside_radius(position));
    
    // Resize indices and values
    indices.resize(number_of_functions_);
    vals.resize(number_of_functions_);
    grad_vals.assign(number_of_functions_,
                     vector<double>(dimension_));
    
    // Get indices of basis functions
    for (int i = 0; i < number_of_functions_; ++i)
    {
        indices[i] = neighbor_functions_[i]->index();
    }
    
    // Get meshless function values
    vector<double> base_values(number_of_functions_);
    vector<vector<double> > grad_base_values(number_of_functions_,
                                             vector<double>(dimension_));
    for (int i = 0; i < number_of_functions_; ++i)
    {
        base_values[i] = neighbor_functions_[i]->value(position);
        grad_base_values[i] = neighbor_functions_[i]->gradient_value(position);
    }
     
    // Get center positions
    vector<vector<double> > center_positions(number_of_functions_);
    for (int i = 0; i < number_of_functions_; ++i)
    {
        center_positions[i] = neighbor_functions_[i]->position();
    }
    
    // Get values
    normalization_->get_gradient_values(position,
                                        center_positions,
                                        base_values,
                                        grad_base_values,
                                        vals,
                                        grad_vals);
}

void Quadratic_MLS_Function::
output(XML_Node output_node) const
{
    output_node.set_attribute("linear_mls", "type");
    function_->output(output_node);
}

void Quadratic_MLS_Function::
check_class_invariants() const
{
    Assert(dimension_ >= 1);
    Assert(number_of_functions_ >= 4); // number of functions must be greater than number of polynomials (2 for linear)
    Assert(function_);
    for (shared_ptr<Meshless_Function> func : neighbor_functions_)
    {
        AssertMsg(func, "basis for mls function not initialized");
    }
}

void Quadratic_MLS_Function::
get_polynomial(vector<double> const &position,
               vector<double> &poly) const
{
    poly.resize(number_of_polynomials_);
    
    switch (dimension_)
    {
    case 1:
    {
        double x = position[0];
        
        poly[0] = 1;
        poly[1] = x;
        poly[2] = x * x;
        break;
    }
    case 2:
    {
        double x = position[0];
        double y = position[1];
        
        poly[0] = 1;
        poly[1] = x;
        poly[2] = y;
        poly[3] = x * y;
        poly[4] = x * x;
        poly[5] = y * y;
        break;
    }
    case 3:
    {
        double x = position[0];
        double y = position[1];
        double z = position[2];
        
        poly[0] = 1;
        poly[1] = x;
        poly[2] = y;
        poly[3] = z;
        poly[4] = x * y;
        poly[5] = x * z;
        poly[6] = y * z;
        poly[7] = x * x;
        poly[8] = y * y;
        poly[9] = z * z;
        break;
    }
    default:
        AssertMsg(false, "dimension not found");
    }
}

void Quadratic_MLS_Function::
get_d_polynomial(int dim,
                 vector<double> const &position,
                 vector<double> &d_poly) const
{
    d_poly.resize(number_of_polynomials_);

    switch (dimension_)
    {
    case 1:
    {
        double x = position[0];
        
        d_poly[0] = 0;
        d_poly[1] = 1;
        d_poly[2] = 2 * x;
        break;
    }
    case 2:
    {
        double x = position[0];
        double y = position[1];
        
        d_poly[0] = 0.;
        switch (dim)
        {
        case 0:
            d_poly[1] = 1;
            d_poly[2] = 0;
            d_poly[3] = y;
            d_poly[4] = 2 * x;
            d_poly[5] = 0;
            break;
        case 1:
            d_poly[1] = 0;
            d_poly[2] = 1;
            d_poly[3] = x;
            d_poly[4] = 0;
            d_poly[5] = 2 * y;
            break;
        default:
            AssertMsg(false, "dimension not found");
        }
        break;
    }
    case 3:
    {
        double x = position[0];
        double y = position[1];
        double z = position[2];
            
        d_poly[0] = 0;
        switch (dim)
        {
        case 0:
            d_poly[1] = 1;
            d_poly[2] = 0;
            d_poly[3] = 0;
            d_poly[4] = y;
            d_poly[5] = z;
            d_poly[6] = 0;
            d_poly[7] = 2 * x;
            d_poly[8] = 0;
            d_poly[9] = 0;
            break;
        case 1:
            d_poly[1] = 0;
            d_poly[2] = 1;
            d_poly[3] = 0;
            d_poly[4] = x;
            d_poly[5] = 0;
            d_poly[6] = z;
            d_poly[7] = 0;
            d_poly[8] = 2 * y;
            d_poly[9] = 0;
            break;
        case 2:
            d_poly[1] = 0;
            d_poly[2] = 0;
            d_poly[3] = 1;
            d_poly[4] = 0;
            d_poly[5] = x;
            d_poly[6] = y;
            d_poly[7] = 0;
            d_poly[8] = 0;
            d_poly[9] = 2 * z;
            break;
        default:
            AssertMsg(false, "dimension not found");
        }
        break;
    }
    default:
        AssertMsg(false, "dimension not found");
    }
}

void Quadratic_MLS_Function::
get_a(vector<double> const &position,
      vector<double> &a) const
{
    a .assign(number_of_polynomials_ * number_of_polynomials_, 0);
    for (shared_ptr<Meshless_Function> func : neighbor_functions_)
    {
        if (func->inside_radius(position))
        {
            vector<double> poly;
            get_polynomial(func->position(),
                           poly);
            double weight = func->value(position);

            for (int j = 0; j < number_of_polynomials_; ++j)
            {
                for (int k = 0; k < number_of_polynomials_; ++k)
                {
                    int l = j + number_of_polynomials_ * k;
                    double polyprod = poly[j] * poly[k];
                    a[l] += weight * polyprod;
                }
            }
        }
    }
}

void Quadratic_MLS_Function::
get_d_a(int dim,
        vector<double> const &position,
        vector<double> &a,
        vector<double> &d_a) const
{
    a.assign(number_of_polynomials_ * number_of_polynomials_, 0);
    d_a.assign(number_of_polynomials_ * number_of_polynomials_, 0);
    
    for (shared_ptr<Meshless_Function> func : neighbor_functions_)
    {
        if (func->inside_radius(position))
        {
            vector<double> poly;
            get_polynomial(func->position(),
                           poly);
            double weight = func->value(position);
            double d_weight = func->d_value(dim,
                                            position);
                      
            for (int j = 0; j < number_of_polynomials_; ++j)
            {
                for (int k = 0; k < number_of_polynomials_; ++k)
                {
                    int l = j + number_of_polynomials_ * k;
                    double polyprod = poly[j] * poly[k];
                    a[l] += weight * polyprod;
                    d_a[l] += d_weight * polyprod;
                }
            }
        }
    }
}

void Quadratic_MLS_Function::
get_b(vector<double> const &position,
      vector<double> &b) const
{
    b.resize(number_of_polynomials_);
    
    vector<double> poly;
    get_polynomial(function_->position(),
                   poly);
    double weight = function_->value(position);
    
    for (int i = 0; i < number_of_polynomials_; ++i)
    {
        b[i] = poly[i] * weight;
    }
}

void Quadratic_MLS_Function::
get_d_b(int dim,
        vector<double> const &position,
        vector<double> &b,
        vector<double> &d_b) const
{
    b.resize(number_of_polynomials_);
    d_b.resize(number_of_polynomials_);
    
    vector<double> poly;
    get_polynomial(function_->position(),
                   poly);
    double weight = function_->value(position);
    double d_weight = function_->d_value(dim,
                                         position);
    
    for (int i = 0; i < number_of_polynomials_; ++i)
    {
        b[i] = poly[i] * weight;
        d_b[i] = poly[i] * d_weight;
    }
}

shared_ptr<Meshless_Normalization> Quadratic_MLS_Function::
normalization() const
{
    return normalization_;
}
