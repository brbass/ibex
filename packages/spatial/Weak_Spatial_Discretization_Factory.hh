#ifndef Weak_Spatial_Discretization_Factory_hh
#define Weak_Spatial_Discretization_Factory_hh

#include <memory>
#include <string>
#include <vector>

#include "Weight_Function.hh"

class Basis_Function;
class Cartesian_Plane;
class Constructive_Solid_Geometry;
class Distance;
class KD_Tree;
class Meshless_Function;
class RBF;
class Weak_Spatial_Discretization;

class Weak_Spatial_Discretization_Factory
{
public:
    
    // Constructor
    Weak_Spatial_Discretization_Factory(std::shared_ptr<Constructive_Solid_Geometry> solid_geometry);
    
    // Get Cartesian grid of points from Cartesian planes
    void get_cartesian_points(std::vector<int> dimensional_points,
                              int &number_of_points,
                              std::vector<std::vector<double> > &points) const;
    
    // Get KD tree from a set of points
    std::shared_ptr<KD_Tree> get_kd_tree(int number_of_points,
                                         std::vector<std::vector<double> > &points) const;
    
    // Get neighbors for each point, given a radius
    // Assumes points for the main set and the neighbor set
    // share center positions and therefore kd_tree.
    void get_neighbors(std::shared_ptr<KD_Tree> kd_tree,
                       int number_of_points,
                       std::vector<double> const &radii,
                       std::vector<double> const &other_radii,
                       std::vector<std::vector<double> > const &positions,
                       std::vector<std::vector<int> > &neighbors) const;
    
    // Get independent meshless functions
    void get_rbf_functions(int number_of_points,
                           std::vector<double> const &radii,
                           std::vector<std::vector<double> > const &points,
                           std::shared_ptr<RBF> rbf,
                           std::shared_ptr<Distance> distance,
                           std::vector<std::shared_ptr<Meshless_Function> > &functions) const;

    // Get MLS functions
    void get_mls_functions(int number_of_points,
                           std::vector<std::shared_ptr<Meshless_Function> > const &functions,
                           std::vector<std::vector<int> > const &neighbors,
                           std::vector<std::shared_ptr<Meshless_Function> > &mls_functions) const;
    
    // Get intersecting boundary surfaces for a given meshless function
    void get_boundary_surfaces(std::shared_ptr<Meshless_Function> function,
                               std::vector<std::shared_ptr<Cartesian_Plane> > &local_boundaries) const;

    // Get basis functions
    void get_basis_functions(int number_of_points,
                             std::vector<std::shared_ptr<Meshless_Function> > const &functions,
                             std::vector<std::shared_ptr<Basis_Function> > &bases) const;

    // Get weight functions
    void get_weight_functions(int number_of_points,
                              Weight_Function::Options options,
                              std::vector<std::vector<int> > const &neighbors,
                              std::vector<std::shared_ptr<Meshless_Function> > const &functions,
                              std::vector<std::shared_ptr<Basis_Function> > const &bases,
                              std::vector<std::shared_ptr<Weight_Function> > &weights) const;

    // Get a spatial discretization with constant radii
    std::shared_ptr<Weak_Spatial_Discretization> get_simple_discretization(int num_dimensional_points,
                                                                           double radius_num_intervals,
                                                                           bool basis_mls,
                                                                           bool weight_mls,
                                                                           std::string basis_type,
                                                                           std::string weight_type,
                                                                           Weight_Function::Options options) const;
    
private:

    std::shared_ptr<Constructive_Solid_Geometry> solid_geometry_;
};

#endif
