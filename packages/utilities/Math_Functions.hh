#ifndef Math_Functions_hh
#define Math_Functions_hh

namespace Math_Functions
{
    // Factorial
    int factorial(int n);
    
    // Returns Legendre polynomial
    double legendre_polynomial(int l,
                               double const &x);

    // Returns associated real Legendre polynomial function
    double legendre_polynomial(int l,
                               int m,
                               double const &x);

    // Returns real spherical harmonic function
    double spherical_harmonic_rec(int l,
                                  int m,
                                  double const &mu,
                                  double const &phi);

    double spherical_harmonic_rec(int l,
                                  int m,
                                  double const &x,
                                  double const &y,
                                  double const &z);

    double spherical_harmonic(int l,
                              int m,
                              double const &x,
                              double const &y,
                              double const &z);
    
    void spherical_to_xyz(double const &mu,
                          double const &phi,
                          double &x,
                          double &y,
                          double &z);

    void xyz_to_spherical(double const &x,
                          double const &y,
                          double const &z,
                          double &mu,
                          double &phi);
} // namespace Math_Functions

#endif
