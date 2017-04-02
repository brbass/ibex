#ifndef Dimensional_Moments_hh
#define Dimensional_Moments_hh

#include <vector>

/*
  Describes the dimensional moments that appear in the SUPG discretization
*/
class Dimensional_Moments
{
public:
    
    Dimensional_Moments(bool supg,
                        int dimension);

    int dimension() const
    {
        return dimension_;
    }
    
    // Dimension + 1
    int number_of_dimensional_moments() const
    {
        return number_of_dimensional_moments_;
    }

    // Number of moments, not repeating like terms
    int number_of_double_dimensional_moments() const
    {
        return number_of_double_dimensional_moments_;
    }

    // Go from d,d' to k:
    // k = indices[d + num_dimensional_moments * d']
    std::vector<int> const &dimensional_indices() const
    {
        return dimensional_indices_;
    }

    // Go from k to d,d':
    // d = subscripts[0 + 2 * k]
    // d' = subscripts[1 + 2 * k]
    std::vector<int> const &dimensional_subscripts() const
    {
        return dimensional_subscripts_;
    }

    // Get coefficients for summation
    // In 3D, {1, tau * direction[0], tau * direction[1], tau * direction[2]}
    std::vector<double> coefficients(double tau,
                                     std::vector<double> const &direction) const;
    std::vector<double> double_coefficients(double tau,
                                            std::vector<double> const &direction) const;
private:

    bool supg_;
    int dimension_;
    int number_of_dimensional_moments_;
    int number_of_double_dimensional_moments_;
    std::vector<int> dimensional_indices_;
    std::vector<int> dimensional_subscripts_;
};

#endif
