#ifndef Constructive_Solid_Geometry_Parser_hh
#define Constructive_Solid_Geometry_Parser_hh

#include <vector>

#include "Parser.hh"
#include "Constructive_Solid_Geometry.hh"

using std::vector;

class Constructive_Solid_Geometry_Parser : public Parser<Constructive_Solid_Geometry>
{
public:

    // Creator
    Constructive_Solid_Geometry_Parser(pugi::xml_node &input_file,
                          vector<shared_ptr<Material> > materials,
                          vector<shared_ptr<Boundary_Source> > boundary_sources);

    // Return pointer to object
    virtual shared_ptr<Constructive_Solid_Geometry> get_ptr() override
    {
        return solid_;
    }

    
private:

    // Get solid geometry
    shared_ptr<Constructive_Solid_Geometry> get_solid(pugi::xml_node &solid_node);
    vector<shared_ptr<Surface> > get_surfaces(pugi::xml_node &surfaces_node,
                                              int dimension);
    vector<shared_ptr<Region> > get_regions(pugi::xml_node &regions_node,
                                            vector<shared_ptr<Surface> > const &surfaces);
    
    shared_ptr<Constructive_Solid_Geometry> solid_;
    vector<shared_ptr<Material> > materials_;
    vector<shared_ptr<Boundary_Source> > boundary_sources_;
};

#endif
