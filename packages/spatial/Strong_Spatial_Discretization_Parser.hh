#ifndef Strong_Spatial_Discretization_Parser_hh
#define Strong_Spatial_Discretization_Parser_hh

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
class Strong_Spatial_Discretization;
class Weak_Spatial_Discretization_Options;
class Weight_Function;
class Weight_Function_Options;
class XML_Node;

/*
  Based on the Weak_Spatial_Discretization_Parser
  Get a Weak_Spatial_Discretization that behaves like a strong spatial discretization
*/
class Strong_Spatial_Discretization_Parser
{
public:

    // Constructor
    Strong_Spatial_Discretization_Parser(std::shared_ptr<Solid_Geometry> solid_geometry,
                                         std::vector<std::shared_ptr<Cartesian_Plane> > const &boundary_surfaces);

    // Parse from XML node
    std::shared_ptr<Strong_Spatial_Discretization> get_strong_discretization(XML_Node input_node) const;
    // std::shared_ptr<Weak_Spatial_Discretization> get_cartesian_discretization(XML_Node input_node) const;
    std::shared_ptr<Strong_Spatial_Discretization> get_points_discretization(XML_Node input_node) const;
    std::shared_ptr<Weight_Function_Options> get_weight_options(XML_Node input_node) const;
    std::shared_ptr<Weak_Spatial_Discretization_Options> get_weak_options(XML_Node input_node) const;
    
private:

    double weight_radii_;
    std::shared_ptr<Solid_Geometry> solid_geometry_;
    std::vector<std::shared_ptr<Cartesian_Plane> > boundary_surfaces_;
};

#endif
