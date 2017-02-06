#ifndef Weak_Spatial_Discretization_Parser_hh
#define Weak_Spatial_Discretization_Parser_hh

class Weak_Spatial_Discretization_Parser
{
public:

    // Constructor
    Weak_Spatial_Discretization_Parser();

    // Parse from XML node
    std::shared_ptr<Weak_Spatial_Discretization> get_weak_discretization(XML_Node input_node);
    std::vector<std::shared_ptr<RBF_Function> > get_rbf_functions(XML_Node input_node,
                                                                  int number_of_points,
                                                                  int dimension,
                                                                  std::string prefix);
    std::vector<std::shared_ptr<Linear_MLS_Function> > get_mls_functions(XML_Node input_node,
                                                                         int number_of_points,
                                                                         int dimension,
                                                                         std::string prefix);
    std::vector<std::shared_ptr<Meshless_Function> > get_meshless_functions(XML_Node input_node,
                                                                            int number_of_points,
                                                                            int dimension,
                                                                            std::string prefix);
    std::vector<std::shared_ptr<Basis_Function> > get_basis_functions(XML_Node input_node,
                                                                      int number_of_points,
                                                                      int dimension);
    std::vector<std::shared_ptr<Weight_Function> > get_weight_functions(XML_Node input_node,
                                                                        int number_of_points,
                                                                        int dimension,
                                                                        std::vector<std::shared_ptr<Basis_Function> > const &basis_functions);
    std::vector<std::shared_ptr<Cartesian_Surface> > get_boundary_surfaces(std::shared_ptr<Meshless_Function>);
    
private:

    std::shared_ptr<Solid_Geometry> solid_geometry_;
    std::vector<std::shared_ptr<Cartesian_Plane> > boundary_surfaces_;
};

#endif
