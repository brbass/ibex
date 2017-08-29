#include "Quadrature_Rule.hh"

#include <cmath>
#include <iostream>
#include <functional>
#include <limits>
#include <memory>
#include <vector>

#include "quadrule.hh"

#include "Check.hh"

namespace Quadrature_Rule
{
    using namespace std;
    
    bool gauss_legendre(int n,
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

        return true;
    }

    bool quadrature_1d(Quadrature_Type quadrature_type,
                       int n,
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
        return true;
    }
    
    bool cartesian_1d(Quadrature_Type quadrature_type,
                      int n,
                      double x1,
                      double x2,
                      vector<double> &ordinates,
                      vector<double> &weights)
    {
        if (x1 > x2)
        {
            cerr << "cartesian_1d: x1 > x2" << endl;
            ordinates.resize(0);
            weights.resize(0);
            return false;
        }
        
        // Get quadrature set
        quadrature_1d(quadrature_type,
                      n,
                      ordinates,
                      weights);
        
        // Scale quadrature set
        double dx = x2 - x1;
        double xt = x2 + x1;
        
        for (int i = 0; i < n; ++i)
        {
            ordinates[i] = 0.5 * (xt + dx * ordinates[i]);
            weights[i] = 0.5 * dx * weights[i];
        }

        return true;
    }

    bool cartesian_2d(Quadrature_Type quadrature_type_x,
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
        if (x1 > x2)
        {
            cerr << "cartesian_2d: x1 > x2" << endl;
            ordinates_x.resize(0);
            ordinates_y.resize(0);
            weights.resize(0);
            return false;
        }
        if (y1 > y2)
        {
            cerr << "cartesian_2d: y1 > y2" << endl;
            ordinates_x.resize(0);
            ordinates_y.resize(0);
            weights.resize(0);
            return false;
        }

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
        
        return true;
    }

    bool cartesian_3d(Quadrature_Type quadrature_type_x,
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
        if (x1 > x2)
        {
            cerr << "cartesian_2d: x1 > x2" << endl;
            ordinates_x.resize(0);
            ordinates_y.resize(0);
            ordinates_z.resize(0);
            weights.resize(0);
            return false;
        }
        if (y1 > y2)
        {
            cerr << "cartesian_2d: y1 > y2" << endl;
            ordinates_x.resize(0);
            ordinates_y.resize(0);
            ordinates_z.resize(0);
            weights.resize(0);
            return false;
        }
        if (z1 > z2)
        {
            cerr << "cartesian_3d: z1 > z2" << endl;
            ordinates_x.resize(0);
            ordinates_y.resize(0);
            ordinates_z.resize(0);
            weights.resize(0);
            return false;
        }
        
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
        
        return true;
    }

    bool cylindrical_2d(Quadrature_Type quadrature_type_r,
                        Quadrature_Type quadrature_type_t,
                        int nr,
                        int nt,
                        double x0,
                        double y0,
                        double r1,
                        double r2,
                        double t1,
                        double t2,
                        vector<double> &ordinates_x,
                        vector<double> &ordinates_y,
                        vector<double> &weights)
    {
        if (r1 > r2)
        {
            cerr << "cylindrical_2d: r1 > r2" << endl;
            ordinates_x.resize(0);
            ordinates_y.resize(0);
            weights.resize(0);
            return false;
        }
        if (t1 > t2)
        {
            cerr << "cylindrical_2d: r1 > r2" << endl;
            ordinates_x.resize(0);
            ordinates_y.resize(0);
            weights.resize(0);
            return false;
        }

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
        
        return true;
    }
    
    bool spherical_3d(Quadrature_Type quadrature_type_r,
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
        if (r1 > r2)
        {
            cerr << "spherical_3d: r1 > r2" << endl;
            ordinates_x.resize(0);
            ordinates_y.resize(0);
            ordinates_z.resize(0);
            weights.resize(0);
            return false;
        }
        if (t1 > t2)
        {
            cerr << "spherical_3d: t1 > t2" << endl;
            ordinates_x.resize(0);
            ordinates_y.resize(0);
            ordinates_z.resize(0);
            weights.resize(0);
            return false;
        }
        if (f1 > f2)
        {
            cerr << "spherical_3d: f1 > f2" << endl;
            ordinates_x.resize(0);
            ordinates_y.resize(0);
            ordinates_z.resize(0);
            weights.resize(0);
            return false;
        }


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
        
        return true;
    }

    bool double_cylindrical_2d(Quadrature_Type quadrature_type_xi,
                               Quadrature_Type quadrature_type_eta,
                               int nxi,
                               int neta,
                               double x1,
                               double y1,
                               double r1,
                               double x2,
                               double y2,
                               double r2,
                               vector<double> &ordinates_x,
                               vector<double> &ordinates_y,
                               vector<double> &weights)
    {
        // Geometric data

        double dx = x2 - x1;
        double dy = y2 - y1;
        double d = sqrt(dx * dx + dy * dy);

        // Find type of quadrature

        if (d >= r1 + r2) // circles do not have any common points
        {
            cerr << "double_cylindrical_2d: no intersections" << endl;
            ordinates_x.resize(0);
            ordinates_y.resize(0);
            weights.resize(0);
            return false;
        }
        else if (d < abs(r1 - r2) + 1e-15) // one circle contained in other
        {
            if (r1 < r2) // first circle contained in second
            {
                return cylindrical_2d(quadrature_type_xi,
                                      quadrature_type_eta,
                                      nxi,
                                      neta,
                                      x1,
                                      y1,
                                      0,
                                      r1,
                                      0,
                                      2 * M_PI,
                                      ordinates_x,
                                      ordinates_y,
                                      weights);
            }
            else // second circle contained in first
            {
                return cylindrical_2d(quadrature_type_xi,
                                      quadrature_type_eta,
                                      nxi,
                                      neta,
                                      x2,
                                      y2,
                                      0,
                                      r2,
                                      0,
                                      2 * M_PI,
                                      ordinates_x,
                                      ordinates_y,
                                      weights);
            }
        }
        else //lens-shaped region
        {
            // Get 1D quadrature sets
            
            vector<double> ord_xi;
            vector<double> ord_eta;
            vector<double> wei_xi;
            vector<double> wei_eta;
            
            quadrature_1d(quadrature_type_xi,
                          nxi,
                          ord_xi,
                          wei_xi);
            quadrature_1d(quadrature_type_eta,
                          neta,
                          ord_eta,
                          wei_eta);
            
            double x_intercept = (d * d + r1 * r1 - r2 * r2) / (2 * d);
            double sqrt_val = -(d - r1 - r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2);
            double y_intercept = sqrt(sqrt_val) / (2 * d);
            
            // Get 2D quadrature set
        
            int n = nxi * neta;
            ordinates_x.resize(n);
            ordinates_y.resize(n);
            weights.resize(n);
            
            if (x_intercept < 0) // intercept past midpoint of first circle
            {
                for (int j = 0; j < neta; ++j)
                {
                    double eta = ord_eta[j];
                    double y_tilde = r1 * eta;
                    double a = (abs(y_tilde) > y_intercept
                                ? -sqrt(r1 * r1 - y_tilde * y_tilde)
                                : d - sqrt(r2 * r2 - y_tilde * y_tilde));
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
                        weights[k] =  0.5 * (b - a) * r1 * wei_xi[i] * wei_eta[j];
                    }
                }
                
            }
            else if (x_intercept > d) // intercept past midpoint of second circle
            {
                for (int j = 0; j < neta; ++j)
                {
                    double eta = ord_eta[j];
                    double y_tilde = r2 * eta;
                    double a = d - sqrt(r2 * r2 - y_tilde * y_tilde);
                    double b = (abs(y_tilde) > y_intercept
                                ? d + sqrt(r2 * r2 - y_tilde * y_tilde)
                                : sqrt(r1 * r1 - y_tilde * y_tilde));
                    
                    for (int i = 0; i < nxi; ++i)
                    {
                        double xi = ord_xi[i];
                        double x_tilde = 0.5 * (b - a) * xi + 0.5 * (a + b);
                        double x = x1 + (dx * x_tilde - dy * y_tilde) / d;
                        double y = y1 + (dy * x_tilde + dx * y_tilde) / d;
                        
                        int k = j + neta * i;
                        
                        ordinates_x[k] = x;
                        ordinates_y[k] = y;
                        weights[k] =  0.5 * (b - a) * r2 * wei_xi[i] * wei_eta[j];
                    }
                }
            }
            else // normal, lens-shaped region
            {
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
            }
            
            return true;
        }
    }

    bool cartesian_bounded_cylindrical_2d(Quadrature_Type quadrature_type_xi,
                                          Quadrature_Type quadrature_type_eta,
                                          int nxi,
                                          int neta,
                                          double x0,
                                          double y0,
                                          double r,
                                          double xminb,
                                          double xmaxb,
                                          double yminb,
                                          double ymaxb,
                                          vector<double> &ordinates_x,
                                          vector<double> &ordinates_y,
                                          vector<double> &weights)
    {
        // Check whether the area is nonzero
        
        double xminc = x0 - r;
        double xmaxc = x0 + r;
        double yminc = y0 - r;
        double ymaxc = y0 + r;

        if (xmaxb < xminc || xminb > xmaxc || ymaxb < yminc || yminb > xmaxc)
        {
            cerr << "cartesian_bounded_cylindrial_2d: no intersections" << endl;
            ordinates_x.resize(0);
            ordinates_y.resize(0);
            weights.resize(0);
            return false;
        }
        if (xminb > xmaxb)
        {
            cerr << "cartesian_bounded_cylindrial_2d: xminb > xmaxb" << endl;
            ordinates_x.resize(0);
            ordinates_y.resize(0);
            weights.resize(0);
            return false;
        }
        if (yminb > ymaxb)
        {
            cerr << "cartesian_bounded_cylindrial_2d: yminb > ymaxb" << endl;
            ordinates_x.resize(0);
            ordinates_y.resize(0);
            weights.resize(0); 
            return false;
        }

        // Check for intersection with boundaries: if none, return cylindrical
        
        if (x0-r > xminb && x0+r < xmaxb && y0-r > yminb && y0+r < ymaxb)
        {
            return cylindrical_2d(quadrature_type_xi,
                                  quadrature_type_eta,
                                  nxi,
                                  neta,
                                  x0,
                                  y0,
                                  0, /*r1*/
                                  r,
                                  0, /*t1*/
                                  2 * M_PI, /*t2*/
                                  ordinates_x,
                                  ordinates_y,
                                  weights);
        }
        
        // Get 1D quadratures
        
        vector<double> ord_xi;
        vector<double> ord_eta;
        vector<double> wei_xi;
        vector<double> wei_eta;
        
        quadrature_1d(quadrature_type_xi,
                      nxi,
                      ord_xi,
                      wei_xi);
        quadrature_1d(quadrature_type_eta,
                      neta,
                      ord_eta,
                      wei_eta);
        
        // Get 2D quadrature
        
        int n = nxi * neta;
        ordinates_x.resize(n);
        ordinates_y.resize(n);
        weights.resize(n);
        
        double ymin = max(yminb, yminc); // yminb > yminc ? yminb : yminc;
        double ymax = min(ymaxb, ymaxc); // ymaxb < ymaxc ? ymaxb : ymaxc;
        
        for (int j = 0; j < neta; ++j)
        {
            double eta = ord_eta[j];
            double y = 0.5 * (ymax * (1 + eta) + ymin * (1 - eta));
            double dy = y -y0;
            double sqrtry = sqrt(r * r - dy * dy);
            double xminc = x0 - sqrtry;
            double xmaxc = x0 + sqrtry;
            double xmin = max(xminb, xminc); // xminb > xminc ? xminb : xminc;
            double xmax = min(xmaxb, xmaxc); // xmaxb < xmaxc ? xmaxb : xmaxc;
            
            for (int i = 0; i < nxi; ++i)
            {
                double xi = ord_xi[i];
                double x = 0.5 * (xmax * (1 + xi) + xmin * (1 - xi));
                
                int k = j + neta * i;
                
                ordinates_x[k] = x;
                ordinates_y[k] = y;
                weights[k] = 0.25 * (xmax - xmin) * (ymax - ymin) * wei_xi[i] * wei_eta[j];
            }
        }
        
        return true;
    }
                                  
    bool cartesian_bounded_double_cylindrical_2d(Quadrature_Type quadrature_type_xi,
                                                 Quadrature_Type quadrature_type_eta,
                                                 int nxi,
                                                 int neta,
                                                 double x1,
                                                 double y1,
                                                 double r1,
                                                 double x2,
                                                 double y2,
                                                 double r2,
                                                 double xminb,
                                                 double xmaxb,
                                                 double yminb,
                                                 double ymaxb,
                                                 vector<double> &ordinates_x,
                                                 vector<double> &ordinates_y,
                                                 vector<double> &weights)
    {
        double x1min = x1 - r1;
        double x1max = x1 + r1;
        double y1min = y1 - r1;
        double y1max = y1 + r1;
        double x2min = x2 - r2;
        double x2max = x2 + r2;
        double y2min = y2 - r2;
        double y2max = y2 + r2;

        double xmin = max(x1min, x2min);
        double xmax = min(x1max, x2max);
        double ymin = max(y1min, y2min);
        double ymax = min(y1max, y2max);

        // Return normal double_cylindrical if no boundaries intersect
        if (xmin > xminb && xmax < xmaxb && ymin > yminb && ymax < ymaxb)
        {
            return double_cylindrical_2d(quadrature_type_xi,
                                         quadrature_type_eta,
                                         nxi,
                                         neta,
                                         x1,
                                         y1,
                                         r1,
                                         x2,
                                         y2,
                                         r2,
                                         ordinates_x,
                                         ordinates_y,
                                         weights);
        }
        
        xmin = max(xminb, xmin);
        xmax = min(xmaxb, xmax);
        ymin = max(yminb, ymin);
        ymax = min(ymaxb, ymax);
        
        return cartesian_2d(quadrature_type_xi,
                            quadrature_type_eta,
                            nxi,
                            neta,
                            xmin,
                            xmax,
                            ymin,
                            ymax,
                            ordinates_x,
                            ordinates_y,
                            weights);
    }
    
    void convert_to_position_1d(vector<double> const &ordinates_x,
                                vector<vector<double> > &ordinates)
    {
        int dimension = 1;
        int number_of_ordinates = ordinates_x.size();
        vector<double> position_template(dimension, 0);
        ordinates.assign(number_of_ordinates, position_template);
        for (int i = 0; i < number_of_ordinates; ++i)
        {
            ordinates[i][0] = ordinates_x[i];
        }
    }
    
    void convert_to_position_2d(vector<double> const &ordinates_x,
                                vector<double> const &ordinates_y,
                                vector<vector<double> > &ordinates)
    {
        int dimension = 2;
        int number_of_ordinates = ordinates_x.size();
        vector<double> position_template(dimension, 0);
        ordinates.assign(number_of_ordinates, position_template);
        for (int i = 0; i < number_of_ordinates; ++i)
        {
            ordinates[i][0] = ordinates_x[i];
            ordinates[i][1] = ordinates_y[i];
        }
    }
    void convert_to_position_3d(vector<double> const &ordinates_x,
                                vector<double> const &ordinates_y,
                                vector<double> const &ordinates_z,
                                vector<vector<double> > &ordinates)
    {
        int dimension = 3;
        int number_of_ordinates = ordinates_x.size();
        vector<double> position_template(dimension, 0);
        ordinates.assign(number_of_ordinates, position_template);
        for (int i = 0; i < number_of_ordinates; ++i)
        {
            ordinates[i][0] = ordinates_x[i];
            ordinates[i][1] = ordinates_y[i];
            ordinates[i][2] = ordinates_z[i];
        }
    }
}
