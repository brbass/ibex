#ifndef Boundary_Source_Parser_hh
#define Boundary_Source_Parser_hh

#include <memory>
#include <vector>

class Angular_Discretization;
class Boundary_Source;
class Energy_Discretization;
class XML_Node;

/*
  Create a Boundary_Source object from XML file
*/
class Boundary_Source_Parser
{
public:

    // Constructor
    Boundary_Source_Parser(std::shared_ptr<Angular_Discretization> angular,
                           std::shared_ptr<Energy_Discretization> energy);
    
    // Parse from XML node
    std::vector<std::shared_ptr<Boundary_Source> > parse_from_xml(XML_Node input_file);
    std::shared_ptr<Boundary_Source> get_vacuum_boundary();
    std::shared_ptr<Boundary_Source> get_reflective_boundary();
    
private:
    
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
};

#endif
