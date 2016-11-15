#include "Quadrature_Rule.hh"

#include <cmath>
#include <iostream>
#include <limits>
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
    
    void cartesian_1d(Quadrature_Type quadrature_type,
                      int n,
                      double x1,
                      double x2,
                      vector<double> &ordinates,
                      vector<double> &weights)
    {
        // Get quadrature set
        switch(quadrature_type)
        {
        case Quadrature_Type::GAUSS_LEGENDRE:
            gauss_legendre(n,
                           ordinates,
                           weights);
            break;
        }
        
        // Scale quadrature set
        double dx = x2 - x1;
        double xt = x2 + x1;
        
        for (int i = 0; i < n; ++i)
        {
            ordinates[i] = 0.5 * (xt + dx * ordinates[i]);
            weights[i] = 0.5 * dx * weights[i];
        }

        return;
    }

    void cartesian_2d(Quadrature_Type quadrature_type_x,
                      Quadrature_Type quadrature_type_y,
                      int nx,
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
        
        cartesian_1d(quadrature_type_x,
                     nx,
                     x1,
                     x2,
                     ord_x,
                     wei_x);
        cartesian_1d(quadrature_type_y,
                     ny,
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

    void cartesian_3d(Quadrature_Type quadrature_type_x,
                      Quadrature_Type quadrature_type_y,
                      Quadrature_Type quadrature_type_z,
                      int nx,
                      int ny,
                      int nz,
                      double x1,
                      double x2,
                      double y1,
                      double y2,
                      double z1,
                      double z2,
                      std::vector<double> &ordinates_x,
                      std::vector<double> &ordinates_y,
                      std::vector<double> &ordinates_z,
                      std::vector<double> &weights)
    {
        // Get 1D quadrature sets
        
        vector<double> ord_x;
        vector<double> ord_y;
        vector<double> ord_z;
        vector<double> wei_x;
        vector<double> wei_y;
        vector<double> wei_z;
        
        cartesian_1d(quadrature_type_x,
                     nx,
                     x1,
                     x2,
                     ord_x,
                     wei_x);
        cartesian_1d(quadrature_type_y,
                     ny,
                     y1,
                     y2,
                     ord_y,
                     wei_y);
        cartesian_1d(quadrature_type_z,
                     nz,
                     z1,
                     z2,
                     ord_z,
                     wei_z);
        
        // Fill 2D quadrature sets
        
        int n = nx * ny * nz;
        ordinates_x.resize(n);
        ordinates_y.resize(n);
        ordinates_z.resize(n);
        weights.resize(n);
        
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int k = 0; k < nz; ++k)
                {
                    int l = k + nz * (j + ny * i);
                
                    ordinates_x[l] = ord_x[i];
                    ordinates_y[l] = ord_y[j];
                    ordinates_z[l] = ord_z[k];
                    weights[l] = wei_x[i] * wei_y[j] * wei_z[k];
                }
            }
        }
        
        return;
    }

    void cylindrical_2d(Quadrature_Type quadrature_type_r,
                        Quadrature_Type quadrature_type_t,
                        int nr,
                        int nt,
                        double x0,
                        double y0,
                        double r1,
                        double r2,
                        double t1,
                        double t2,
                        std::vector<double> &ordinates_x,
                        std::vector<double> &ordinates_y,
                        std::vector<double> &weights)
    {
        // Get 1D quadrature sets
        
        vector<double> ord_r;
        vector<double> ord_t;
        vector<double> wei_r;
        vector<double> wei_t;
        
        cartesian_1d(quadrature_type_r,
                     nr,
                     r1,
                     r2,
                     ord_r,
                     wei_r);
        cartesian_1d(quadrature_type_t,
                     nt,
                     t1,
                     t2,
                     ord_t,
                     wei_t);

        // Fill 2D quadrature sets
        
        int n = nr * nt;
        ordinates_x.resize(n);
        ordinates_y.resize(n);
        weights.resize(n);
        
        for (int i = 0; i < nr; ++i)
        {
            double r = ord_r[i];
            
            for (int j = 0; j < nt; ++j)
            {
                int k = j + nt * i;
                double t = ord_t[j];
                
                ordinates_x[k] = x0 + r * cos(t);
                ordinates_y[k] = y0 + r * sin(t);
                weights[k] = wei_r[i] * wei_t[j] * r;
            }
        }
        
        return;
    }
    
    void spherical_3d(Quadrature_Type quadrature_type_r,
                      Quadrature_Type quadrature_type_t,
                      Quadrature_Type quadrature_type_f,
                      int nr,
                      int nt,
                      int nf,
                      double x0,
                      double y0,
                      double z0,
                      double r1,
                      double r2,
                      double t1,
                      double t2,
                      double f1,
                      double f2,
                      std::vector<double> &ordinates_x,
                      std::vector<double> &ordinates_y,
                      std::vector<double> &ordinates_z,
                      std::vector<double> &weights)
    {
        // Get 1D quadrature sets
        
        vector<double> ord_r;
        vector<double> ord_t;
        vector<double> ord_f;
        vector<double> wei_r;
        vector<double> wei_t;
        vector<double> wei_f;
        
        cartesian_1d(quadrature_type_r,
                     nr,
                     r1,
                     r2,
                     ord_r,
                     wei_r);
        cartesian_1d(quadrature_type_t,
                     nt,
                     t1,
                     t2,
                     ord_t,
                     wei_t);
        cartesian_1d(quadrature_type_f,
                     nf,
                     f1,
                     f2,
                     ord_f,
                     wei_f);
        
        // Fill 3D quadrature sets
        
        int n = nr * nt * nf;
        ordinates_x.resize(n);
        ordinates_y.resize(n);
        ordinates_z.resize(n);
        weights.resize(n);
        
        for (int i = 0; i < nr; ++i)
        {
            double r = ord_r[i];
            
            for (int j = 0; j < nt; ++j)
            {
                double t = ord_t[j];
                
                for (int k = 0; k < nf; ++k)
                {
                    int l = k + nf * (j + nt * i);
                    double f = ord_f[k];
                    
                    ordinates_x[l] = x0 + r * cos(t) * sin(f);
                    ordinates_y[l] = y0 + r * sin(t) * sin(f);
                    ordinates_z[l] = z0 + r * cos(f);
                    weights[l] = wei_r[i] * wei_t[j] * wei_f[k] * r * r * sin(f);
                }
            }
        }
        
        return;
    }

    bool lens_2d(Quadrature_Type quadrature_type_xi,
                 Quadrature_Type quadrature_type_eta,
                 int nxi,
                 int neta,
                 double x1,
                 double y1,
                 double x2,
                 double y2,
                 double r1,
                 double r2,
                 vector<double> &ordinates_x,
                 vector<double> &ordinates_y,
                 vector<double> &weights)
    {
        // Get 1D quadrature sets
        
        vector<double> ord_xi;
        vector<double> ord_eta;
        vector<double> wei_xi;
        vector<double> wei_eta;

        gauss_legendre(nxi,
                       ord_xi,
                       wei_xi);
        gauss_legendre(neta,
                       ord_eta,
                       wei_eta);
        
        // Geometric data
        
        double dx = x2 - x1;
        double dy = y2 - y1;
        double d = sqrt(dx * dx + dy * dy);
        
        if (d == 0)
        {
            cerr << "lens_2d: distance between centers is zero" << endl;
            return false;
        }
        
        double x_intercept = (d * d + r1 * r1 - r2 * r2) / (2 * d);
        if (x_intercept < 0)
        {
            cerr << "lens_2d: shape is not a lens" << endl;
        }
        double sqrt_val = 2 * d * d * (r1 * r1 + r2 * r2) - pow(r1 * r1 - r2 * r2, 2) - pow(d, 4);
        if (sqrt_val <= 0)
        {
            cerr << "lens_2d: no intersection" << endl;
            return false;
        }
        double y_intercept = sqrt(sqrt_val) / (2 * d);

        // Get 2D quadrature set

        int n = nxi * neta;
        ordinates_x.resize(n);
        ordinates_y.resize(n);
        weights.resize(n);
        
        for (int j = 0; j < neta; ++j)
        {
            double eta = ord_eta[j];
            double y_tilde = y_intercept * eta;
            double a = d - sqrt(r2 * r2 - y_tilde * y_tilde);
            double b = sqrt(r1 * r1 - y_tilde * y_tilde);

            for (int i = 0; i < nxi; ++i)
            {
                double xi = ord_xi[i];
                double x_tilde = 0.5 * (b - a) * xi + 0.5 * (a + b);
                double x = x1 + (dx * x_tilde - dy * y_tilde) / d;
                double y = y1 + (dy * x_tilde + dx * y_tilde) / d;
                
                int k = j + neta * i;
                
                ordinates_x[k] = x;
                ordinates_y[k] = y;
                weights[k] =  0.5 * (b - a) * y_intercept * wei_xi[i] * wei_eta[j];
            }
        }
        
        return true;
    }
        
}
