#ifndef Quadrature_Rule_hh
#define Quadrature_Rule_hh

#include <vector>

/*
  Returns quadratures in x (1D), x,y (2D) and x,y,z (3D)
*/
namespace Quadrature_Rule
{
    /*
      Quadrature types
    */

    enum class Quadrature_Type
    {
        GAUSS_LEGENDRE
    };

    /*
      Get vectors of Gauss_Legendre ordinates and weights
    */
    
    void gauss_legendre(int n, std::vector<double> &ordinates, std::vector<double> &weights);
    
    /*
      1D Cartesian quadrature
    */

    void cartesian_1d(Quadrature_Type quadrature_type,
                      int n,
                      double x1,
                      double x2,
                      std::vector<double> &ordinates,
                      std::vector<double> &weights);
    
    /*
      2D Cartesian quadrature
    */

    void cartesian_2d(Quadrature_Type quadrature_type_x,
                      Quadrature_Type quadrature_type_y,
                      int nx, 
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
                      std::vector<double> &weights);

    /*
      2D cylindrical quadrature
    */

    void cylindrical_2d(Quadrature_Type quadrature_type_r,
                        Quadrature_Type quadrature_type_t,
                        int nr,
                        int nt,
                        double x0,
                        double y0,
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
                      std::vector<double> &ordinates_r,
                      std::vector<double> &ordinates_t,
                      std::vector<double> &ordinates_f,
                      std::vector<double> &weights);
}
#endif


