#include "Dimensional_Moments.hh"

using std::vector;

Dimensional_Moments::
Dimensional_Moments(bool supg,
                    int dimension):
    supg_(supg),
    dimension_(dimension)
{
    // Get packed symmetric storage
    // Simpler to hardcode than use formulas
    if (supg)
    {
        number_of_dimensional_moments_ = dimension_ + 1;
        switch (dimension)
        {
        case 1:
            number_of_double_dimensional_moments_ = 3;
            dimensional_indices_
                = {0, 1,
                   1, 2};
            dimensional_subscripts_
                = {0, 0,
                   0, 1,
                   1, 1};
            break;
        case 2:
            number_of_double_dimensional_moments_ = 6;
            dimensional_indices_
                = {0, 1, 2,
                   1, 3, 4,
                   2, 4, 5};
            dimensional_subscripts_
                = {0, 0,
                   0, 1,
                   0, 2,
                   1, 1,
                   1, 2,
                   2, 2};
            break;
        case 3:
            number_of_double_dimensional_moments_ = 10;
            dimensional_indices_
                = {0, 1, 2, 3,
                   1, 4, 5, 6,
                   2, 5, 7, 8,
                   3, 6, 8, 9};
            dimensional_subscripts_
                = {0, 0,
                   0, 1,
                   0, 2,
                   0, 3,
                   1, 1,
                   1, 2,
                   1, 3,
                   2, 2,
                   2, 3,
                   3, 3};
            break;
        }
    }
    else
    {
        number_of_dimensional_moments_ = 1;
        number_of_double_dimensional_moments_ = 1;
        dimensional_indices_ = {0};
        dimensional_subscripts_ = {0, 0};
    }
}

vector<double> Dimensional_Moments::
coefficients(double tau,
             vector<double> const &direction) const
{
    vector<double> vals(number_of_dimensional_moments_);
    
    vals[0] = 1;

    for (int d = 1; d < number_of_dimensional_moments_; ++d)
    {
        vals[d] = tau * direction[d - 1];
    }

    return vals;
}

vector<double> Dimensional_Moments::
double_coefficients(double tau,
                    vector<double> const &direction) const
{
    vector<double> coeffs = coefficients(tau,
                                         direction);
    vector<double> vals(number_of_double_dimensional_moments_);
    for (int d = 0; d < number_of_double_dimensional_moments_; ++d)
    {
        int d1 = dimensional_subscripts_[0 + 2 * d];
        int d2 = dimensional_subscripts_[1 + 2 * d];
        
        vals[d] = coeffs[d1] * coeffs[d2];
    }

    return vals;
}
