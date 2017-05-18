#include "Linear_MLS_External.hh"

#include "Linear_Algebra.hh"
#include "Matrix_Functions.hh"
#include "Vector_Functions.hh"

namespace la = Linear_Algebra;
namespace mf = Matrix_Functions;
namespace vf = Vector_Functions;
using std::vector;

Linear_MLS_External::
Linear_MLS_External(int dimension):
    dimension_(dimension),
    number_of_polynomials_(dimension + 1)
{
}

void Linear_MLS_External::
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

void Linear_MLS_External::
get_grad_polynomial(vector<double> const &position,
                    vector<vector<double> > &grad_poly) const
{
    grad_poly.assign(dimension_, vector<double>(dimension_ + 1, 0));

    for (int d = 0; d < dimension_; ++d)
    {
        grad_poly[d][d + 1] = 1.;
    }
}

void Linear_MLS_External::
get_values(vector<double> const &position,
           vector<vector<double> > const &center_positions,
           vector<double> const &base_values,
           vector<double> &values) const
{
    // Get size information
    int number_of_functions = base_values.size();
    values.resize(number_of_functions);
    
    // Get polynomial value at evaluation position
    vector<double> p_position(number_of_polynomials_);
    get_polynomial(position,
                   p_position);

    // Get polynomial values at centers
    vector<vector<double> > p_center(number_of_functions, vector<double>(number_of_polynomials_));
    for (int i = 0; i < number_of_functions; ++i)
    {
        get_polynomial(center_positions[i],
                       p_center[i]);
    }
    
    // Get A matrix
    vector<double> a_mat(number_of_polynomials_ * number_of_polynomials_);
    for (int i = 0; i < number_of_functions; ++i)
    {
        double const weight = base_values[i];
        vector<double> const &poly = p_center[i];
        for (int j = 0; j < number_of_polynomials_; ++j)
        {
            for (int k = 0; k < number_of_polynomials_; ++k)
            {
                int l = j + number_of_polynomials_ * k;
                a_mat[l] += weight * poly[j] * poly[k];
            }
        }
    }
    
    // Get inverse of A
    vector<double> a_inv_mat(number_of_polynomials_ * number_of_polynomials_);
    la::direct_inverse(a_mat,
                       a_inv_mat);
    
    // Get values of basis function
    for (int i = 0; i < number_of_functions; ++i)
    {
        // Get B vector
        vector<double> b_vec(number_of_polynomials_);
        {
            double const weight = base_values[i];
            vector<double> const &poly = p_center[i];
            for (int j = 0; j < number_of_polynomials_; ++j)
            {
                b_vec[j] = weight * poly[j];
            }
        }

        // Get value of function
        values[i] = vf::dot(p_position,
                            mf::square_matrix_vector_product(number_of_polynomials_,
                                                             a_inv_mat,
                                                             b_vec));
    }
    
}
           

void Linear_MLS_External::
get_gradient_values(vector<double> const &position,
                    vector<vector<double> > const &center_positions,
                    vector<double> const &base_values,
                    vector<vector<double> > const &grad_base_values,
                    vector<double> &values,
                    vector<vector<double> > &gradient_values) const
{
    // Get size information
    int number_of_functions = base_values.size();
    values.resize(number_of_functions);
    gradient_values.assign(number_of_functions, vector<double>(dimension_, 0));

    // Get polynomial value at evaluation position
    vector<double> p_position(number_of_polynomials_);
    get_polynomial(position,
                   p_position);
    vector<vector<double> > grad_p_position(dimension_,
                                            vector<double>(number_of_polynomials_));
    get_grad_polynomial(position,
                        grad_p_position);

    
    // Get polynomial values at centers
    vector<vector<double> > p_center(number_of_functions, vector<double>(number_of_polynomials_));
    for (int i = 0; i < number_of_functions; ++i)
    {
        get_polynomial(center_positions[i],
                       p_center[i]);
    }
    
    // Get A matrix
    vector<double> a_mat(number_of_polynomials_ * number_of_polynomials_);
    vector<vector<double> > grad_a_mat(dimension_,
                                       vector<double>(number_of_polynomials_ * number_of_polynomials_));
    for (int i = 0; i < number_of_functions; ++i)
    {
        double const weight = base_values[i];
        vector<double> const &poly = p_center[i];
        for (int j = 0; j < number_of_polynomials_; ++j)
        {
            for (int k = 0; k < number_of_polynomials_; ++k)
            {
                int l = j + number_of_polynomials_ * k;
                a_mat[l] += weight * poly[j] * poly[k];
            }
        }

        for (int d = 0; d < dimension_; ++d)
        {
            double const d_weight = grad_base_values[i][d];
            for (int j = 0; j < number_of_polynomials_; ++j)
            {
                for (int k = 0; k < number_of_polynomials_; ++k)
                {
                    int l = j + number_of_polynomials_ * k;
                    grad_a_mat[d][l] += d_weight * poly[j] * poly[k];
                }
            }
        }
    }
    
    // Get inverse of A
    vector<double> a_inv_mat(number_of_polynomials_ * number_of_polynomials_);
    la::direct_inverse(a_mat,
                       a_inv_mat);
    vector<vector<double> > grad_a_inv_mat(dimension_,
                                           vector<double>(number_of_polynomials_ * number_of_polynomials_));
    for (int d = 0; d < dimension_; ++d)
    {
        vector<double> &mat = grad_a_inv_mat[d];
        mat = mf::square_matrix_matrix_product(number_of_polynomials_,
                                               grad_a_mat[d],
                                               a_inv_mat);
        mat = mf::square_matrix_matrix_product(number_of_polynomials_,
                                               a_inv_mat,
                                               mat);
        mat = vf::multiply(mat,
                           -1.);
    }
    
    // Get values of basis function
    for (int i = 0; i < number_of_functions; ++i)
    {
        // Get B vector
        vector<double> b_vec(number_of_polynomials_);
        double const weight = base_values[i];
        vector<double> const &poly = p_center[i];
        for (int j = 0; j < number_of_polynomials_; ++j)
        {
            b_vec[j] = weight * poly[j];
        }
        
        // Get value of function
        values[i] = vf::dot(p_position,
                            mf::square_matrix_vector_product(number_of_polynomials_,
                                                             a_inv_mat,
                                                             b_vec));

        for (int d = 0; d < dimension_; ++d)
        {
            // Get d_B vector
            vector<double> d_b_vec(number_of_polynomials_);
            double const d_weight = grad_base_values[i][d];

            for (int j = 0; j < number_of_polynomials_; ++j)
            {
                d_b_vec[j] = d_weight * poly[j];
            }
            
            // Get gradient values of function

            double t1 = vf::dot(grad_p_position[d],
                                mf::square_matrix_vector_product(number_of_polynomials_,
                                                                 a_inv_mat,
                                                                 b_vec));
            double t2 = vf::dot(p_position,
                                mf::square_matrix_vector_product(number_of_polynomials_,
                                                                 grad_a_inv_mat[d],
                                                                 b_vec));
            double t3 = vf::dot(p_position,
                                mf::square_matrix_vector_product(number_of_polynomials_,
                                                                 a_inv_mat,
                                                                 d_b_vec));
            gradient_values[i][d] = t1 + t2 + t3;
        }
    }
}    
