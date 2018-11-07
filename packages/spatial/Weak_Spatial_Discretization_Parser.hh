#ifndef Weak_Spatial_Discretization_Parser_hh
#define Weak_Spatial_Discretization_Parser_hh

#include <memory>
#include <vector>

class Basis_Function;
class Cartesian_Plane;
class Dimensional_Moments;
class Distance;
class Linear_MLS_Function;
class Meshless_Function;
class RBF;
class RBF_Function;
class Solid_Geometry;
class Weak_Spatial_Discretization;
class Weak_Spatial_Discretization_Options;
class Weight_Function;
class Weight_Function_Options;
class XML_Node;

class Weak_Spatial_Discretization_Parser
{
public:

    // Constructor
    Weak_Spatial_Discretization_Parser(std::shared_ptr<Solid_Geometry> solid_geometry,
                                       std::vector<std::shared_ptr<Cartesian_Plane> > const &boundary_surfaces);

    // Parse from XML node
    std::shared_ptr<Weak_Spatial_Discretization> get_weak_discretization(XML_Node input_node) const;
    std::shared_ptr<Weak_Spatial_Discretization> get_full_discretization(XML_Node input_node) const;
    std::shared_ptr<Weak_Spatial_Discretization> get_cartesian_discretization(XML_Node input_node) const;
    std::shared_ptr<Weak_Spatial_Discretization> get_galerkin_points_discretization(XML_Node input_node) const;
    std::shared_ptr<Weak_Spatial_Discretization> get_legendre_discretization(XML_Node input_node) const;
    std::shared_ptr<Weight_Function_Options> get_weight_options(XML_Node input_node) const;
    std::shared_ptr<Weak_Spatial_Discretization_Options> get_weak_options(XML_Node input_node) const;
    std::vector<std::shared_ptr<Meshless_Function> > get_rbf_functions(XML_Node input_node,
                                                                       int number_of_points,
                                                                       int dimension,
                                                                       std::string prefix) const;
    std::vector<std::shared_ptr<Meshless_Function> > get_mls_functions(XML_Node input_node,
                                                                       int order,
                                                                       int number_of_points,
                                                                       int dimension,
                                                                       std::string prefix) const;
    std::vector<std::shared_ptr<Meshless_Function> > get_meshless_functions(XML_Node input_node,
                                                                            int number_of_points,
                                                                            int dimension,
                                                                            std::string prefix) const;
    std::vector<std::shared_ptr<Basis_Function> > get_basis_functions(XML_Node input_node,
                                                                      int number_of_points,
                                                                      int dimension) const;
    std::vector<std::shared_ptr<Weight_Function> > get_weight_functions(XML_Node input_node,
                                                                        int number_of_points,
                                                                        int dimension,
                                                                        std::shared_ptr<Weak_Spatial_Discretization_Options> options,
                                                                        std::shared_ptr<Dimensional_Moments> dimensional_moments,
                                                                        std::vector<std::shared_ptr<Basis_Function> > const &basis_functions) const;
    std::vector<std::shared_ptr<Cartesian_Plane> > get_boundary_surfaces(std::shared_ptr<Meshless_Function> function) const;
    
private:

    std::shared_ptr<Solid_Geometry> solid_geometry_;
    std::vector<std::shared_ptr<Cartesian_Plane> > boundary_surfaces_;
};

#endif
