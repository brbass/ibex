#ifndef Weak_Spatial_Discretization_Factory_hh
#define Weak_Spatial_Discretization_Factory_hh

#include <memory>
#include <string>
#include <vector>

#include "Meshless_Function_Factory.hh"

class Basis_Function;
class Cartesian_Plane;
class Dimensional_Moments;
class KD_Tree;
class Meshless_Function;
class Solid_Geometry;
class Weak_Spatial_Discretization;
class Weak_Spatial_Discretization_Options;
class Weight_Function;
class Weight_Function_Options;

class Weak_Spatial_Discretization_Factory
{
public:
    
    // Constructor
    Weak_Spatial_Discretization_Factory(std::shared_ptr<Solid_Geometry> solid_geometry,
                                        std::vector<std::shared_ptr<Cartesian_Plane> > const &boundary_surfaces);
    
    // Get basis functions
    void get_basis_functions(int number_of_points,
                             std::vector<std::shared_ptr<Meshless_Function> > const &functions,
                             std::vector<std::shared_ptr<Basis_Function> > &bases) const;

    // Get weight functions
    void get_weight_functions(int number_of_points,
                              std::shared_ptr<Weight_Function_Options> weight_options,
                              std::shared_ptr<Weak_Spatial_Discretization_Options> weak_options,
                              std::shared_ptr<Dimensional_Moments> dimensional_moments,
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
                                                                           std::shared_ptr<Weight_Function_Options> weight_options,
                                                                           std::shared_ptr<Weak_Spatial_Discretization_Options> weak_options) const;
    
private:

    std::shared_ptr<Solid_Geometry> solid_geometry_;
    std::vector<std::shared_ptr<Cartesian_Plane> > boundary_surfaces_;
    Meshless_Function_Factory meshless_factory_;
};

#endif
