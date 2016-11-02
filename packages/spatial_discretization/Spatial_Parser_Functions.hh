#ifndef Spatial_Parser_Functions_hh
#define Spatial_Parser_Functions_hh

#include <memory>
#include <vector>

class Boundary_Source;
class Cartesian_Overlay;
class Distance;
class Material;
class Constructive_Solid_Geometry;

using std::shared_ptr;
using std::vector;

namespace Spatial_Parser_Functions
{
    void get_point(int dimension,
                   double bounding_radius,
                   vector<double> const &bounding_origin,
                   vector<double> &point);
    
    void get_ray(int dimension,
                 double bounding_radius,
                 vector<double> const &bounding_origin,
                 vector<double> &origin,
                 vector<double> &direction);
    
    void get_random_points(shared_ptr<Constructive_Solid_Geometry> solid_geometry,
                           shared_ptr<Distance> distance_metric,
                           int dimension,
                           int max_attempts,
                           double min_distance_boundary,
                           double min_distance_internal,
                           double bounding_radius,
                           vector<double> const &bounding_origin,
                           int &number_of_points,
                           int &number_of_boundary_points,
                           int &number_of_internal_points,
                           vector<int> &surfaces,
                           vector<int> &regions,
                           vector<int> &boundary_points,
                           vector<int> &internal_points,
                           vector<vector<double> > &positions,
                           vector<vector<double> > &boundary_normal,
                           vector<shared_ptr<Material> > &material,
                           vector<shared_ptr<Boundary_Source> > &boundary_source,
                           shared_ptr<Cartesian_Overlay> cartesian_overlay);
    
    void get_neighbor_information(int index,
                                  int dimension,
                                  int number_of_points,
                                  int number_of_groups,
                                  int number_of_neighbors,
                                  int number_to_average,
                                  shared_ptr<Distance> const distance_metric,
                                  vector<vector<double> > const &positions,
                                  vector<int> &neighbor_indices,
                                  vector<double> &mean_distance);
    
    void get_shape_parameter(int number_of_groups,
                             double shape_multiplier,
                             vector<double> &mean_distance,
                             vector<double> &shape_parameter);

    void get_optical_min_distance(int number_of_groups,
                                  double min_cartesian_distance,
                                  vector<double> const &sigma_t,
                                  vector<double> &min_distance);                         
} // namespace Spatial_Parser_Functions

#endif
