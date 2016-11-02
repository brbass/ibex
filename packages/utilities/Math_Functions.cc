#include "Math_Functions.hh"

#include "Check.hh"

#include <cmath>
#include <cstdlib>

namespace Math_Functions
{
    using namespace std;
    
    int factorial(int n)
    {
        int val = 1;
        
        for (int i = n; i >= 1; --i)
        {
            val *= i; 
        }
        
        return val;

        // return (n == 0) ? 1 : n * factorial(n - 1);
    }
    
    int factorial2(int n)
    {
        int val = 1;
        
        for (int i = n; i >=1; i -= 2)
        {
            val *= i;
        }

        return val;
    }

    double legendre_polynomial(int l,
                               double const &x)
    {
        if (l == 0)
        {
            return 1;
        }
        
        double plm2 = 0;
        double plm1 = 1;
        double pl = x;
    
        for (int i = 2; i <= l; ++i)
        {
            double j = static_cast<double>(i);
        
            plm2 = plm1;
            plm1 = pl;
            pl = ((2 * j - 1) * x * plm1 - (j - 1) * plm2) / j;
        }
        
        return pl;
    }
    
    double legendre_polynomial(int l,
                               int m,
                               double const &x)
    {
        if (l == m)
        {
            return factorial2(2 * l - 1) * pow(1. - x * x, l / 2.);
        }
        else if (l == m + 1)
        {
            return factorial2(2 * l - 1) * x * pow(1. - x * x, (l - 1.) / 2.);
        }

        double p2 = 0;
        double p1 = factorial2(2 * m - 1) * pow(1. - x * x, m / 2.);
        double p = x * (2 * m + 1) * p1;
        
        for (int i = m + 2; i <= l; ++i)
        {
            p2 = p1;
            p1 = p;
            
            p =  (x * (2 * i - 1.) * p1 - (i + m - 1.) * p2) / (i - m);
        }
        
        return p;
    }

    double spherical_harmonic_rec(int l,
                                  int m,
                                  double const &mu,
                                  double const &phi)
    {
        double val = (m == 0) ? 1 : 0;
        
        val = sqrt((2 - val) * factorial(l - abs(m)) / factorial(l + abs(m)));
        
        double t = (m >= 0) ? cos(m * phi) : sin(abs(m) * phi);
        
        return val * legendre_polynomial(l, abs(m), mu) * t;
    }

    double spherical_harmonic_rec(int l,
                                  int m,
                                  double const &x,
                                  double const &y,
                                  double const &z)
    {
        double mu;
        double phi;
        
        xyz_to_spherical(x, y, z, mu, phi);
        
        return spherical_harmonic_rec(l, m, mu, phi);
    }

    double spherical_harmonic(int l,
                              int m,
                              double const &x,
                              double const &y,
                              double const &z)
    {
        switch(l)
        {
        case 0:
            return 1;
        case 1:
            switch(m)
            {
            case -1:
                return z;
            case 0:
                return x;
            case 1:
                return y;
            default:
                AssertMsg(false, "Could not find moment");
                break;
            }
        case 2:
            switch(m)
            {
            case -2:
                return sqrt(3.) * y * z;
            case -1:
                return sqrt(3.) * x * z;
            case 0:
                return 0.5 * (3. * x * x - 1.);
            case 1:
                return sqrt(3.) * x * y;
            case 2:
                return 0.5 * sqrt(3.) * (y * y - z * z);
            default:
                AssertMsg(false, "Could not find moment");
                break;
            }
        case 3:
            switch(m)
            {
            case -3:
                return sqrt(0.625) * z * (3. * y * y - z * z);
            case -2:
                return sqrt(15.) * x * y * z;
            case -1:
                return sqrt(0.375) * z * (5. * x * x - 1);
            case 0:
                return 0.5 * x * (5. * x * x - 3.);
            case 1:
                return sqrt(0.375) * y * (5. * x * x - 1.);
            case 2:
                return sqrt(3.75) * x * (y * y - z * z);
            case 3:
                return sqrt(0.625) * y * (y * y - 3 * z * z);
            default:
                AssertMsg(false, "Could not find moment");
                break;
            }
        default:
            return spherical_harmonic_rec(l, m, x, y, z);
        }
    }
    
    // See file sphere_lebedev_rule.cpp for original conversion code
    void xyz_to_spherical(double const &x,
                          double const &y,
                          double const &z,
                          double &mu,
                          double &phi)
    {
        double val = sqrt(y * y + z * z);
        
        if (val > 0)
        {
            phi = acos(y / val);
        }
        else
        {
            phi = acos(y);
        }
        
        if (z < 0)
        {
            phi = -phi;
        }
        
        mu = x;
    }
    
    void spherical_to_xyz(double const &mu,
                          double const &phi,
                          double &x,
                          double &y,
                          double &z)
    {
        double val = sqrt(1 - mu * mu);
        
        x = mu;
        y = cos(phi) * val;
        z = sin(phi) * val;
    }

} // namespace Math_Functions
