#include "Angular_Discretization.hh"

#include <cmath>
#include <limits>
#include <vector>

#include "Check.hh"
#include "Math_Functions.hh"

using namespace std;

Angular_Discretization::
Angular_Discretization(int dimension, 
                       int number_of_scattering_moments,
                       int number_of_ordinates):
    dimension_(dimension),
    number_of_scattering_moments_(number_of_scattering_moments),
    number_of_ordinates_(number_of_ordinates)
{
    initialize_moment_data();
}

void Angular_Discretization::
initialize_moment_data()
{
    switch(dimension_)
    {
    case 1:
    {
        angular_normalization_ = 2;
        
        int m = 0;
        
        for (int l = 0; l < number_of_scattering_moments_; ++l)
        {
            l_indices_.push_back(l);
            m_indices_.push_back(m);
        }
        
        number_of_moments_ = number_of_scattering_moments_;
        
        break;
    }
    case 2:
    {
        angular_normalization_ = 2 * M_PI;
        
        int sum = 0;
        for (int l = 0; l < number_of_scattering_moments_; ++l)
        {
            for (int m = 0; m <= l; ++m)
            {
                l_indices_.push_back(l);
                m_indices_.push_back(m);
                
                sum += 1;
            }
        }
        
        number_of_moments_ = sum;
        
        break;
    }
    case 3:
    {
        angular_normalization_ = 4 * M_PI;
        
        int sum = 0;
        for (int l = 0; l < number_of_scattering_moments_; ++l)
        {
            for (int m = -l; m <= l; ++m)
            {
                l_indices_.push_back(l);
                m_indices_.push_back(m);
                
                sum += 1;
            }
        }
        
        number_of_moments_ = sum;
        
        break;
    }
    default:
        AssertMsg(false, "Dimension not found");
        break;
    }
}

double Angular_Discretization::
moment(int mom,
       int ord) const
{
    switch(dimension_)
    {
    case 1:
    {
        double mu = ordinates()[ord];
        
        return Math_Functions::legendre_polynomial(mom, mu);
    }
    case 2:
    {
        int l = l_indices_[mom];
        int m = m_indices_[mom];
        
        int o_x = 0 + dimension_ * ord;
        int o_y = 1 + dimension_ * ord;
        double x = ordinates()[o_x];
        double y = ordinates()[o_y];
        double z = sqrt(1 - x * x - y * y);

        return Math_Functions::spherical_harmonic(l, m, x, y, z);
    }
    case 3:
    {
        int l = l_indices_[mom];
        int m = m_indices_[mom];
        
        int o_x = 0 + dimension_ * ord;
        int o_y = 1 + dimension_ * ord;
        int o_z = 2 + dimension_ * ord;
        double x = ordinates()[o_x];
        double y = ordinates()[o_y];
        double z = ordinates()[o_z];
        
        return Math_Functions::spherical_harmonic(l, m, x, y, z);
    }
    default:
        AssertMsg(false, "dimension not found");
        return 0;
    }
}

void Angular_Discretization::
moment_to_discrete(vector<double> &data) const
{
    // Check size of input vector
    Check(data.size() == number_of_moments_);
    
    // Make copy of input vector
    vector<double> old_data(data);

    // Resize input vector
    data.resize(number_of_ordinates_);

    // Perform moment to discrete operation
    for (int o = 0; o < number_of_ordinates_; ++o)
    {
        double sum = 0;
        
        for (int m = 0; m < number_of_moments_; ++m)
        {
            int const l = scattering_indices()[m];
            double const p = moment(m, o);
            
            sum += (2 * static_cast<double>(l) + 1) / angular_normalization_ * p * old_data[m];
        }

        data[o] = sum;
    }
}

void Angular_Discretization::
discrete_to_moment(vector<double> &data) const
{
    // Check size of input vector
    Check(data.size() == number_of_ordinates_);
    
    // Make copy of input vector
    vector<double> old_data(data);

    // Resize input vector
    data.resize(number_of_moments_);

    // Perform moment to discrete operation
    for (int m = 0; m < number_of_moments_; ++m)
    {
        double sum = 0;
        for (int o = 0; o < number_of_ordinates_; ++o)
        {
            double const p = moment(m, o);
            
            sum += weights()[o] * p * old_data[o];
        }

        data[m] = sum;
    }
}

void Angular_Discretization::
manufactured_coefficients(vector<int> &size,
                          vector<vector<int> > &indices,
                          vector<vector<double> > &coefficients) const
{
    // Initialize full coefficient matrix to zero
    vector<vector<double> > coeff_matrix(number_of_moments_ * number_of_moments_, vector<double>(dimension_, 0));

    // Get moments in all directions
    vector<double> moments(number_of_ordinates_ * number_of_moments_);
    for (int o = 0; o < number_of_ordinates_; ++o)
    {
        for (int m = 0; m < number_of_moments_; ++m)
        {
            int k = m + number_of_moments_ * o;

            moments[k] = moment(m, o);
        }
    }
    
    // Perform integration without normalization term
    for (int o = 0; o < number_of_ordinates_; ++o)
    {
        vector<double> const dir = direction(o);
        double const wei = weights()[o];
        
        for (int m1 = 0; m1 < number_of_moments_; ++m1)
        {
            int const k1 = m1 + number_of_moments_ * o;
            
            for (int m2 = 0; m2 <= m1; ++m2)
            {
                int const k2 = m2 + number_of_moments_ * o;
                int const km = m1 + number_of_moments_ * m2;

                for (int d = 0; d < dimension_; ++d)
                {
                    coeff_matrix[km][d] += wei * dir[d] * moments[k1] * moments[k2];
                }
            }
        }
    }

    // Fill in lower triangle of matrix
    for (int m1 = 0; m1 < number_of_moments_; ++m1)
    {
        for (int m2 = m1 + 1; m2 < number_of_moments_; ++m2)
        {
            int klo = m1 + number_of_moments_ * m2;
            int kup = m2 + number_of_moments_ * m1;
                
            for (int d = 0; d < dimension_; ++d)
            {
                coeff_matrix[klo][d] = coeff_matrix[kup][d];
            }
        }
    }

    // Add in normalization term
    for (int m1 = 0; m1 < number_of_moments_; ++m1)
    {
        int const l = scattering_indices()[m1];
        double const mult = (2 * static_cast<double>(l) + 1) / angular_normalization_;
        for (int m2 = 0; m2 < number_of_moments_; ++m2)
        {
            int k = m1 + number_of_moments_ * m2;
            for (int d = 0; d < dimension_; ++d)
            {
                coeff_matrix[k][d] *= mult;
            }
        }
    }
    
    // Get list of nonzero coefficients
    double tolerance = 1e-12;
    size.resize(number_of_moments_);
    indices.resize(number_of_moments_);
    coefficients.resize(number_of_moments_);
    for (int m2 = 0; m2 < number_of_moments_; ++m2)
    {
        vector<int> moment_ind;
        vector<double> moment_coeffs;
        for (int m1 = 0; m1 < number_of_moments_; ++m1)
        {
            int k = m1 + number_of_moments_ * m2;

            bool include = false;
            for (int d = 0; d < dimension_; ++d)
            {
                if (abs(coeff_matrix[k][d]) > 1e-12)
                {
                    include = true;
                    break;
                }
            }

            if (include)
            {
                moment_ind.push_back(m1);
                for (int d = 0; d < dimension_; ++d)
                {
                    moment_coeffs.push_back(coeff_matrix[k][d]);
                }
            }
        }
        
        size[m2] = moment_ind.size();
        indices[m2] = moment_ind;
        coefficients[m2] = moment_coeffs;
    }
}
