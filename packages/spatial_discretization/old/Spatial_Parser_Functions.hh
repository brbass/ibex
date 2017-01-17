#ifndef Spatial_Parser_Functions_hh
#define Spatial_Parser_Functions_hh

#include <memory>
#include <vector>

class Boundary_Source;
class Cartesian_Overlay;
class Distance;
class Material;
class Constructive_Solid_Geometry;

namespace Spatial_Parser_Functions
{
    void get_point(int dimension,
                   double bounding_radius,
                   std::vector<double> const &bounding_origin,
                   std::vector<double> &point);
    
    void get_ray(int dimension,
                 double bounding_radius,
                 std::vector<double> const &bounding_origin,
                 std::vector<double> &origin,
                 std::vector<double> &direction);
    
    void get_random_points(std::shared_ptr<Constructive_Solid_Geometry> solid_geometry,
                           std::shared_ptr<Distance> distance_metric,
                           int dimension,
                           int max_attempts,
                           double min_distance_boundary,
                           double min_distance_internal,
                           double bounding_radius,
                           std::vector<double> const &bounding_origin,
                           int &number_of_points,
                           int &number_of_boundary_points,
                           int &number_of_internal_points,
                           std::vector<int> &surfaces,
                           std::vector<int> &regions,
                           std::vector<int> &boundary_points,
                           std::vector<int> &internal_points,
                           std::vector<std::vector<double> > &positions,
                           std::vector<std::vector<double> > &boundary_normal,
                           std::vector<std::shared_ptr<Material> > &material,
                           std::vector<std::shared_ptr<Boundary_Source> > &boundary_source,
                           std::shared_ptr<Cartesian_Overlay> cartesian_overlay);
    
    void get_neighbor_information(int index,
                                  int dimension,
                                  int number_of_points,
                                  int number_of_groups,
                                  int number_of_neighbors,
                                  int number_to_average,
                                  std::shared_ptr<Distance> const distance_metric,
                                  std::vector<std::vector<double> > const &positions,
                                  std::vector<int> &neighbor_indices,
                                  std::vector<double> &mean_distance);
    
    void get_shape_parameter(int number_of_groups,
                             double shape_multiplier,
                             std::vector<double> &mean_distance,
                             std::vector<double> &shape_parameter);

    void get_optical_min_distance(int number_of_groups,
                                  double min_cartesian_distance,
                                  std::vector<double> const &sigma_t,
                                  std::vector<double> &min_distance);                         
} // namespace Spatial_Parser_Functions

#endif
