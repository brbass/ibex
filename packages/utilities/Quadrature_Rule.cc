#include "Quadrature_Rule.hh"

#include <iostream>
#include <memory>
#include <vector>

#include "quadrule.hh"

#include "Check.hh"

namespace Quadrature_Rule
{
    using namespace std;
    
    void gauss_legendre(int n,
                        vector<double> &ordinates,
                        vector<double> &weights)
    {
        if (n < 1)
        {
            AssertMsg(false, "quadrature must have n >= 1");
        }

        ordinates.resize(n);
        weights.resize(n);
    
        if (n <= 33
            || (63 <= n && n <= 65)
            || (127 <= n && n <= 129)
            || (255 <= n && n <= 257))
        {
            quadrule::legendre_set(n, &ordinates[0], &weights[0]);
        }
        else
        {
            quadrule::legendre_dr_compute(n, &ordinates[0], &weights[0]);
        }

        return;
    }
    
    void gl_cartesian_1d(int n,
                         double x1,
                         double x2,
                         vector<double> &ordinates,
                         vector<double> &weights)
    {
        // Get quadrature set
        gauss_legendre(n,
                       ordinates,
                       weights);
        
        // Scale quadrature set
        double dx = x2 - x1;
        double xt = x2 + x1;
        
        for (int i = 0; i < order_; ++i)
        {
            ordinates[i] = 0.5 * (xt + dx * ordinates[i]);
            weights[i] = 0.5 * dx * weights[i];
        }

        return;
    }

    void gl_cartesian_2d(int nx,
                         int ny,
                         double x1,
                         double x2,
                         double y1,
                         double y2,
                         vector<double> &ordinates_x,
                         vector<double> &ordinates_y,
                         vector<double> &weights)
    {
        // Get 1D quadrature sets
        
        vector<double> ord_x;
        vector<double> ord_y;
        vector<double> wei_x;
        vector<double> wei_y;
        
        cartesian_1d(nx,
                     x1,
                     x2,
                     ord_x,
                     wei_x);
        cartesian_1d(ny,
                     y1,
                     y2,
                     ord_y,
                     wei_y);
        
        // Fill 2D quadrature sets
        
        int n = nx * ny;
        ordinates_x.resize(n);
        ordinates_y.resize(n);
        weights.resize(n);
        
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                int k = j + ny * i;
                
                ordinates_x[k] = ord_x[i];
                ordinates_y[k] = ord_y[j];
                weights[k] = wei_x[i] * wei_y[j];
            }
        }
        
        return;
    }
}
