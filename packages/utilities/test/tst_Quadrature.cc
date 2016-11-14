#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "Quadrature_Rule.hh"

using namespace std;

namespace qr = Quadrature_Rule;

double compact_gaussian(double epsilon,
                        double s,
                        double smax)
{
    if (s <= smax)
    {
        double e2 = epsilon * epsilon;
        double s2 = s * s;
        double sm2 = smax * smax;
        
        return (exp(-e2 * s2) - exp(-e2 * sm2)) / (1 - exp(-e2 * sm2));
    }
    else
    {
        return 0.;
    }
}

double double_compact_gaussian(double ep1,
                               double ep2,
                               double smax1,
                               double smax2,
                               double x,
                               double y,
                               double dist)
{
    double s1 = sqrt(x * x + y * y);
    double x2 = x - dist;
    double s2 = sqrt(x2 * x2 + y * y);

    return (compact_gaussian(ep1,
                             s1,
                             smax1)
            * compact_gaussian(ep2,
                               s2,
                               smax2));
}

int test_gaussian_2d(int order)
{
    int checksum = 0;
    
    double epsilon = 4;
    double smax = 4 / epsilon;
    double x0 = 0;
    double y0 = 0;
    double tmin = 0;
    double tmax = 2 * M_PI;
    double rmin = 0;
    double rmax = smax;
    double xmin = -smax + x0;
    double xmax = smax + x0;
    double ymin = -smax + y0;
    double ymax = smax + y0;
    qr::Quadrature_Type quad_type = qr::Quadrature_Type::GAUSS_LEGENDRE;
    int num_ordinates = order * order;

    double analytic = M_PI * (1 / (epsilon * epsilon) - smax * smax / (exp(epsilon * epsilon * smax * smax) - 1));

    vector<pair<double, string> > integrals;
    
    // Cartesian 2D
    double cart2d = 0;
    {
        vector<double> ordinates_x;
        vector<double> ordinates_y;
        vector<double> weights;

        qr::cartesian_2d(quad_type,
                         quad_type,
                         order,
                         order,
                         xmin,
                         xmax,
                         ymin,
                         ymax,
                         ordinates_x,
                         ordinates_y,
                         weights);
        
        for (int i = 0; i < num_ordinates; ++i)
        {
            double s = sqrt(ordinates_x[i] * ordinates_x[i] + ordinates_y[i] * ordinates_y[i]);
            
            cart2d += weights[i] * compact_gaussian(epsilon,
                                                    s,
                                                    smax);
        }
    }
    integrals.emplace_back(cart2d, "Cart_2D");
    
    // Cylindrical 2D
    double cyl2d = 0;
    {
        vector<double> ordinates_x;
        vector<double> ordinates_y;
        vector<double> weights;

        qr::cylindrical_2d(quad_type,
                           quad_type,
                           order,
                           order,
                           x0,
                           y0,
                           rmin,
                           rmax,
                           tmin,
                           tmax,
                           ordinates_x,
                           ordinates_y,
                           weights);
        
        for (int i = 0; i < num_ordinates; ++i)
        {
            double s = sqrt(ordinates_x[i] * ordinates_x[i] + ordinates_y[i] * ordinates_y[i]);
            
            cyl2d += weights[i] * compact_gaussian(epsilon,
                                                    s,
                                                    smax);
        }
    }
    integrals.emplace_back(cyl2d, "Cyl_2D");

    int w = 16;
    for (int i = 0; i < integrals.size(); ++i)
    {
        cout << setw(w) << "Gauss_2D";
        cout << setw(w) << integrals[i].second;
        cout << setw(w) << order;
        cout << setw(w) << (analytic - integrals[i].first) / analytic;
        cout << endl;
    }
    
    return checksum;
}

int test_double_gaussian_2d(int order)
{
    int checksum = 0;
    
    double ep1 = 4;
    double ep2 = ep1 / 2;
    double smax1 = 4 / ep1;
    double smax2 = 4 / ep2;
    double dist = 2.5;
    double interx = (dist * dist + smax1 * smax1 - smax2 * smax2) / (2 * dist);
    double intery1 = -sqrt(2 * dist * dist * (smax1 * smax1 + smax2 * smax2) - pow(smax1 * smax1 - smax2 * smax2, 2) - dist * dist * dist * dist) / (2 * dist);
    double intery2 = -intery1;
    double xinter1 = 1.0;
    double xinter2 = 0.5;
    qr::Quadrature_Type quad_type = qr::Quadrature_Type::GAUSS_LEGENDRE;
    int num_ordinates = order * order;

    double numeric = 1.1016707860233166491e-10;
    
    vector<pair<double, string> > integrals;
    
    // Cartesian 2D over weight function
    double cart2d = 0;
    {
        vector<double> ordinates_x;
        vector<double> ordinates_y;
        vector<double> weights;

        qr::cartesian_2d(quad_type,
                         quad_type,
                         order,
                         order,
                         -smax1,
                         smax1,
                         -smax1,
                         smax1,
                         ordinates_x,
                         ordinates_y,
                         weights);
        
        for (int i = 0; i < num_ordinates; ++i)
        {
            double s1 = sqrt(ordinates_x[i] * ordinates_x[i] + ordinates_y[i] * ordinates_y[i]);
            
            cart2d += (weights[i]
                       * double_compact_gaussian(ep1,
                                                 ep2,
                                                 smax1,
                                                 smax2,
                                                 ordinates_x[i],
                                                 ordinates_y[i],
                                                 dist));
        }
    }
    integrals.emplace_back(cart2d, "Cart_2D");
    
    // Cylindrical 2D over weight function
    double cyl2d = 0;
    {
        vector<double> ordinates_x;
        vector<double> ordinates_y;
        vector<double> weights;

        qr::cylindrical_2d(quad_type,
                           quad_type,
                           order,
                           order,
                           interx,
                           0,
                           0,
                           intery2,
                           0,
                           2 * M_PI,
                           ordinates_x,
                           ordinates_y,
                           weights);
        
        for (int i = 0; i < num_ordinates; ++i)
        {
            cyl2d += (weights[i]
                      * double_compact_gaussian(ep1,
                                                ep2,
                                                smax1,
                                                smax2,
                                                ordinates_x[i],
                                                ordinates_y[i],
                                                dist));
        }
    }
    integrals.emplace_back(cyl2d, "Cyl_2D");

    // Cartesian 2D over lens
    double cart2dlens = 0;
    {
        vector<double> ordinates_x;
        vector<double> ordinates_y;
        vector<double> weights;
        
        qr::cartesian_2d(quad_type,
                         quad_type,
                         order,
                         order,
                         xinter2,
                         xinter1,
                         intery1,
                         intery2,
                         ordinates_x,
                         ordinates_y,
                         weights);
        
        for (int i = 0; i < num_ordinates; ++i)
        {
            double s1 = sqrt(ordinates_x[i] * ordinates_x[i] + ordinates_y[i] * ordinates_y[i]);
            
            cart2dlens += (weights[i]
                           * double_compact_gaussian(ep1,
                                                     ep2,
                                                     smax1,
                                                     smax2,
                                                     ordinates_x[i],
                                                     ordinates_y[i],
                                                     dist));
        }
    }
    integrals.emplace_back(cart2dlens, "Cart_2D_lens");
    
    // Cylindrical 2D over lens
    double cyl2dlens = 0;
    {
        vector<double> ordinates_x;
        vector<double> ordinates_y;
        vector<double> weights;

        qr::cylindrical_2d(quad_type,
                           quad_type,
                           order,
                           order,
                           0,
                           0,
                           0,
                           smax1,
                           0,
                           2 * M_PI,
                           ordinates_x,
                           ordinates_y,
                           weights);
        
        for (int i = 0; i < num_ordinates; ++i)
        {
            cyl2dlens += (weights[i]
                          * double_compact_gaussian(ep1,
                                                    ep2,
                                                    smax1,
                                                    smax2,
                                                    ordinates_x[i],
                                                    ordinates_y[i],
                                                    dist));
        }
    }
    integrals.emplace_back(cyl2dlens, "Cyl_2D_lens");
    
    int w = 16;
    for (int i = 0; i < integrals.size(); ++i)
    {
        cout << setw(w) << "double_Gauss_2D";
        cout << setw(w) << integrals[i].second;
        cout << setw(w) << order;
        cout << setw(w) << (numeric - integrals[i].first) / numeric;
        cout << endl;
    }
    
    return checksum;
}

int main()
{
    int checksum = 0;

    int w = 16;
    cout << setw(w) << "integrand";
    cout << setw(w) << "method";
    cout << setw(w) << "order";
    cout << setw(w) << "error";
    cout << endl;
    for (int i = 8; i <= 64; i = i + 8)
    {
        checksum += test_gaussian_2d(i);
    }
    for (int i = 4; i <= 1024; i *= 2)
    {
        checksum += test_double_gaussian_2d(i);
    }
    
    return checksum;
}
