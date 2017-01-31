#ifndef Weak_Spatial_Discretization_Parser_hh
#define Weak_Spatial_Discretization_Parser_hh

class Weak_Spatial_Discretization_Parser
{
public:

    // Constructor
    Weak_Spatial_Discretization_Parser();

    // Parse from XML node
    std::shared_ptr<Weak_Spatial_Discretization> parse_from_xml(XML_Node input_node);
};

#endif
