#ifndef Quadrature_Rule_hh
#define Quadrature_Rule_hh

#include <vector>

namespace Quadrature_Rule
{
    /*
      Get vectors of Gauss_Legendre ordinates and weights
    */
    
    void gauss_legendre(int n, std::vector<double> &ordinates, std::vector<double> &weights);
    
    /*
      1D Cartesian quadrature
    */

    void gl_cartesian_1d(int n,
                         double x1,
                         double x2,
                         std::vector<double> &ordinates,
                         std::vector<double> &weights);
    
    /*
      2D Cartesian quadrature
    */

    void gl_cartesian_2d(int nx, 
                         int ny,
                         double x1,
                         double x2,
                         double y1,
                         double y2,
                         std::vector<double> &ordinates_x,
                         std::vector<double> &ordinates_y,
                         std::vector<double> &weights);
    
    /*
      3D Cartesian quadrature
    */

    void gl_cartesian_3d(int nx,
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
                         std::vector<double> &weights);

    /*
      2D cylindrical quadrature
    */

    void gl_cylindrical_2d(int nr,
                           int nt,
                           double r1,
                           double r1,
                           double t1,
                           double t2,
                           std::vector<double> &ordinates_r,
                           std::vector<double> &ordinates_t,
                           std::vector<double> &weights);
    /*
      3D spherical quadrature
    */
    
}
#endif


