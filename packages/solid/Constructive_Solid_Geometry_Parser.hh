#ifndef Constructive_Solid_Geometry_Parser_hh
#define Constructive_Solid_Geometry_Parser_hh

#include <memory>
#include <vector>

class Boundary_Source;
class Cartesian_Plane;
class Constructive_Solid_Geometry;
class Material;
class Surface;
class Region;
class XML_Node;

class Constructive_Solid_Geometry_Parser
{
public:

    // Creator
    Constructive_Solid_Geometry_Parser(std::vector<std::shared_ptr<Material> > materials,
                                       std::vector<std::shared_ptr<Boundary_Source> > boundary_sources);

    // Get solid geometry
    std::shared_ptr<Constructive_Solid_Geometry> parse_from_xml(XML_Node input_node);
    std::vector<std::shared_ptr<Surface> > get_surfaces(XML_Node input_node,
                                                        int dimension);
    std::vector<std::shared_ptr<Region> > get_regions(XML_Node input_node,
                                                      std::vector<std::shared_ptr<Surface> > const &surfaces);
    
private:
    
    std::vector<std::shared_ptr<Material> > materials_;
    std::vector<std::shared_ptr<Boundary_Source> > boundary_sources_;
};

#endif
