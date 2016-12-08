#ifndef Material_Parser_hh
#define Material_Parser_hh

#include <memory>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;
class Material;
class XML_Node;

/*
  Create a Material object from XML file
*/
class Material_Parser
{
public:

    // Constructor
    Material_Parser(std::shared_ptr<Angular_Discretization> angular,
                    std::shared_ptr<Energy_Discretization> energy);
    
    // Return Material object
    std::vector<std::shared_ptr<Material> > parse_from_xml(XML_Node input_file);
    
private:
    
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
};

#endif
