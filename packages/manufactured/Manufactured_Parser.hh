#ifndef Manufactured_Parser_hh
#define Manufactured_Parser_hh

#include <memory>
#include <vector>

class Angular_Discretization;
class Cartesian_Plane;
class Energy_Discretization;
class Manufactured_Solution;
class Solid_Geometry;
class XML_Node;

class Manufactured_Parser
{
public:
    
    // Creator
    Manufactured_Parser(std::shared_ptr<Angular_Discretization> angular,
                        std::shared_ptr<Energy_Discretization> energy);
    
    // Get solid geometry
    void parse_from_xml(XML_Node input_node,
                        std::shared_ptr<Manufactured_Solution> &solution,
                        std::shared_ptr<Solid_Geometry> &solid,
                        std::vector<std::shared_ptr<Cartesian_Plane> > &boundary_surfaces);
    
private:
    
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
};

#endif
