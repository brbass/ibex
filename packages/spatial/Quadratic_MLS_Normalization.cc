#include "Quadratic_MLS_Normalization.hh"

#include "Dense_Solver_Factory.hh"
#include "Matrix_Functions.hh"
#include "Vector_Functions.hh"

namespace mf = Matrix_Functions;
namespace vf = Vector_Functions;
using std::vector;

Quadratic_MLS_Normalization::
Quadratic_MLS_Normalization(int dimension):
    dimension_(dimension),
    number_of_polynomials_(dimension + 1)
{
    Dense_Solver_Factory factory;
    solver_ = factory.get_solver(number_of_polynomials_);
}

void Quadratic_MLS_Normalization::
get_polynomial(vector<double> const &position,
               vector<double> &poly) const
{
    poly.resize(number_of_polynomials_);
    
    switch (dimension_)
    {
    case 1:
    {
        double x = position[1];
        
        poly[0] = 1;
        poly[1] = x;
        poly[2] = x * x;
        break;
    }
    case 2:
    {
        double x = position[1];
        double y = position[2];
        
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
        double x = position[1];
        double y = position[2];
        double z = position[3];
        
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

void Quadratic_MLS_Normalization::
get_d_polynomial(int dim,
                 vector<double> const &position,
                 vector<double> &d_poly) const
{
    d_poly.resize(number_of_polynomials_);

    switch (dimension_)
    {
    case 1:
    {
        d_poly[0] = 0;
        d_poly[1] = 1;
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

void Quadratic_MLS_Normalization::
get_grad_polynomial(vector<double> const &position,
                    vector<vector<double> > &grad_poly) const
{
    grad_poly.assign(dimension_, vector<double>(number_of_polynomials_, 0));

    for (int d = 0; d < dimension_; ++d)
    {
        get_d_polynomial(d,
                         position,
                         grad_poly[d]);
    }
}

void Quadratic_MLS_Normalization::
get_values(vector<double> const &position,
           vector<vector<double> > const &center_positions,
           vector<double> const &base_values,
           vector<double> &values) const
{
    // Get size information
    int number_of_functions = base_values.size();
    Check(position.size() == dimension_);
    Check(center_positions.size() == number_of_functions);

    // Resize but don't assign in case base_values and values are the same
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
    solver_->inverse(a_mat,
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
           

void Quadratic_MLS_Normalization::
get_gradient_values(vector<double> const &position,
                    vector<vector<double> > const &center_positions,
                    vector<double> const &base_values,
                    vector<vector<double> > const &grad_base_values,
                    vector<double> &values,
                    vector<vector<double> > &gradient_values) const
{
    // Get size information
    int number_of_functions = base_values.size();
    Check(position.size() == dimension_);
    Check(center_positions.size() == number_of_functions);
    Check(grad_base_values.size() == number_of_functions);

    // Resize but don't assign in case values and base values are the same
    values.resize(number_of_functions);
    gradient_values.resize(number_of_functions);
    for (int i = 0; i < number_of_functions; ++i)
    {
        gradient_values[i].resize(dimension_);
    }

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
    solver_->inverse(a_mat,
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
