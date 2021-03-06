#ifndef Meshless_Function_Factory_hh
#define Meshless_Function_Factory_hh

#include <memory>
#include <string>
#include <vector>

class Cartesian_Plane;
class Distance;
class KD_Tree;
class Meshless_Function;
class RBF;

class Meshless_Function_Factory
{
public:

    Meshless_Function_Factory();

    // Get boundary limits from Cartesian planes
    void get_boundary_limits(int dimension,
                             std::vector<std::shared_ptr<Cartesian_Plane> > const &boundaries,
                             std::vector<std::vector<double> > &limits) const;
    
    // Get Cartesian grid of points from limits
    void get_cartesian_points(int dimension,
                              std::vector<int> dimensional_points,
                              std::vector<std::vector<double> > const &limits,
                              int &number_of_points,
                              std::vector<std::vector<double> > &points) const;
    
    // Find "number_of_neighbors" nearest points
    // Radius for each point is the distance to the furthest
    // neighbor times the multiplier
    void get_radii_nearest(std::shared_ptr<KD_Tree> kd_tree,
                           int dimension,
                           int number_of_points,
                           int number_of_neighbors,
                           double radius_multiplier,
                           std::vector<double> &radii) const;

    // Find "number_of_neighbors" nearest points
    // Ensure that radius of each of the neighbors is no
    // larger than the distance to the current point
    void get_radii_coverage(std::shared_ptr<KD_Tree> kd_tree,
                            int dimension,
                            int number_of_points,
                            int number_of_neighbors,
                            double radius_multiplier,
                            std::vector<double> &radii) const;
    
    // Get neighbors for each circular region, given a radius
    // Assumes points for the main set and the neighbor set
    // share center positions and therefore kd_tree.
    void get_neighbors(std::shared_ptr<KD_Tree> kd_tree,
                       bool global_rbf,
                       int dimension,
                       int number_of_points,
                       std::vector<double> const &radii,
                       std::vector<double> const &other_radii,
                       std::vector<std::vector<double> > const &positions,
                       std::vector<std::vector<int> > &neighbors,
                       std::vector<std::vector<double> > &distances) const;

    // Get functions that intersect each point, given a radius
    void get_point_neighbors(std::shared_ptr<KD_Tree> kd_tree,
                             int dimension,
                             int number_of_points,
                             std::vector<double> const &radii,
                             std::vector<std::vector<double> > const &positions,
                             std::vector<std::vector<int> > &neighbors,
                             std::vector<std::vector<double> > &distances) const;
    
    // Get independent meshless functions
    void get_rbf_functions(int number_of_points,
                           std::vector<double> const &radii,
                           std::vector<std::vector<double> > const &points,
                           std::shared_ptr<RBF> rbf,
                           std::shared_ptr<Distance> distance,
                           std::vector<std::shared_ptr<Meshless_Function> > &functions) const;
    
    // Get MLS functions, given meshless functions and neighbors
    // Input and output vector of functions cannot be the same
    void get_mls_functions(int order,
                           int number_of_points,
                           std::vector<std::shared_ptr<Meshless_Function> > const &functions,
                           std::vector<std::vector<int> > const &neighbors,
                           std::vector<std::shared_ptr<Meshless_Function> > &mls_functions) const;

    // Get Legendre tensor product functions, given the order in each dimension
    void get_legendre_functions(int dimension,
                                std::vector<int> const &order,
                                std::vector<std::vector<double> > const &limits,
                                int &number_of_points,
                                std::vector<std::shared_ptr<Meshless_Function> > &functions) const;

    
    // Get MLS functions based on RBF functions
    // Creates and destroys KD tree, so it's advisable not to use this
    // where a KD tree is already being created (for weak discretization)
    void get_cartesian_mls_functions(int order,
                                     int dimension,
                                     double radius_num_intervals,
                                     std::vector<int> dimensional_points,
                                     std::vector<std::vector<double> > const &limits,
                                     std::string rbf_type,
                                     std::vector<std::shared_ptr<Meshless_Function> > &functions) const;

    // Given radii, make certain no two points are too close together
    bool check_point_conditioning(int const number_of_points,
                                  std::vector<double> const &radii,
                                  std::vector<std::vector<int> > const &neighbors,
                                  std::vector<std::vector<double> > const &squared_distances) const;
    
    // Get intersecting boundary surfaces for a given meshless function
    void get_boundary_surfaces(std::shared_ptr<Meshless_Function> function,
                               std::vector<std::shared_ptr<Cartesian_Plane> > const &boundaries,
                               std::vector<std::shared_ptr<Cartesian_Plane> > &local_boundaries) const;
    
private:
};

#endif
