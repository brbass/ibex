#ifndef Weak_Spatial_Discretization_Factory_hh
#define Weak_Spatial_Discretization_Factory_hh

#include <memory>
#include <string>
#include <vector>

#include "Meshless_Function_Factory.hh"
#include "Weight_Function.hh"

class Basis_Function;
class Cartesian_Plane;
class Constructive_Solid_Geometry;
class KD_Tree;
class Meshless_Function;
class Weak_Spatial_Discretization;

class Weak_Spatial_Discretization_Factory
{
public:
    
    // Constructor
    Weak_Spatial_Discretization_Factory(std::shared_ptr<Constructive_Solid_Geometry> solid_geometry);
    
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
    Meshless_Function_Factory meshless_factory_;
};

#endif
