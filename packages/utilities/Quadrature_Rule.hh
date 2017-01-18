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
    bool gauss_legendre(int n, 
                        std::vector<double> &ordinates,
                        std::vector<double> &weights);

    /*
      Get specified 1D quadrature
    */
    bool quadrature_1d(Quadrature_Type quadrature_type,
                       int n,
                       std::vector<double> &ordinates,
                       std::vector<double> &weights);
    
    /*
      1D Cartesian product quadrature
      Integral from x1 to x2
    */
    bool cartesian_1d(Quadrature_Type quadrature_type,
                      int n,
                      double x1,
                      double x2,
                      std::vector<double> &ordinates,
                      std::vector<double> &weights);
    
    /*
      2D Cartesian product quadrature
      Double integral from x1 to x2 and y1 to y2
    */
    bool cartesian_2d(Quadrature_Type quadrature_type_x,
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
      3D Cartesian product quadrature
      Triple integral from x1 to x2, y1 to y2 and z1 to z2
    */
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
                      std::vector<double> &weights);

    /*
      2D cylindrical quadrature
      Double integral from r1 to r2 and t1 to t2 for
      circle centered at x0, y0
      Quadrature points are chosen using cylindrical coordinates
    */
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
                        std::vector<double> &ordinates_x,
                        std::vector<double> &ordinates_y,
                        std::vector<double> &weights);
    
    /*
      3D spherical quadrature
      Triple integral from r1 to r2, t1 to t2 and f1 to f2 for
      sphere centered at x0, y0, z0
      Quadrature points are chosen using spherical coordinates
    */
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
                      std::vector<double> &weights);

    /*
      Quadrature over lens or other circle-circle regions
      Double integral for the intersection of the circles centered at
      x1, y1 and x2, y2 with radii r1 and r2
      Returns:
      - Cylindrical quadrature if one circle is contained in other
      - Lens-based quadrature if the two circles intersect
    */
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
                               std::vector<double> &ordinates_x,
                               std::vector<double> &ordinates_y,
                               std::vector<double> &weights);

    /*
      Cylindrical quadrature with Cartesian bounds
      Double integral for circle centered at x0, y0 with radius r
      with Cartesian boundaries at xminb, xmaxb, yminb and ymaxb
      Returns:
      - Cylindrical quadrature if boundaries don't intersect circle
      - Cartesian quadrature over bounding box if Cartesian boundaries apply
    */
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
                                          std::vector<double> &ordinates_x,
                                          std::vector<double> &ordinates_y,
                                          std::vector<double> &weights);

    /*
      Quadrature over lens with Cartesian bounds
      Double integral over intersection of circles at x1, y1 and x2, y2
      with radii r1 and r2 with Cartesian boundaries at xminb, xmaxb, yminb
      and ymaxb
      Finds bounding boxes for each of the two circles and then finds the
      intersection of these bounding boxes and the Cartesian boundaries
      Returns double cylindrical quadrature if boundaries don't intersect
    */
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
                                                 std::vector<double> &ordinates_x,
                                                 std::vector<double> &ordinates_y,
                                                 std::vector<double> &weights);

    /*
      Convert quadrature from vector<double> for each set of ordinates
      to position-based ordinates of form vector<vector<double> >
    */
    void convert_to_position_1d(std::vector<double> const &ordinates_x,
                                std::vector<std::vector<double> > &ordinates);
    void convert_to_position_2d(std::vector<double> const &ordinates_x,
                                std::vector<double> const &ordinates_y,
                                std::vector<std::vector<double> > &ordinates);
    void convert_to_position_3d(std::vector<double> const &ordinates_x,
                                std::vector<double> const &ordinates_y,
                                std::vector<double> const &ordinates_z,
                                std::vector<std::vector<double> > &ordinates);
}

#endif


