#ifndef Spatial_Discretization_Parser_hh
#define Spatial_Discretization_Parser_hh

#include "Parser.hh"
#include "Spatial_Discretization.hh"

class Boundary_Source;
class Material;
class RBF_Discretization;
class Solid_Geometry;

/*
  Create a Spatial_Discretization object from XML file
*/
class Spatial_Discretization_Parser : public Parser<Spatial_Discretization>
{
public:
    
    // Creator
    Spatial_Discretization_Parser(pugi::xml_node &input_file,
                                  vector<shared_ptr<Material> > const &materials,
                                  vector<shared_ptr<Boundary_Source> > const &boundary_sources);
    
    // Return pointer to object
    virtual shared_ptr<Spatial_Discretization> get_ptr() override
    {
        return spatial_;
    }
    
    // Parse a radial basis function mesh
    shared_ptr<RBF_Discretization> get_rbf(pugi::xml_node &spatial);
    
private:
    
    shared_ptr<Spatial_Discretization> spatial_;
    vector<shared_ptr<Material> > materials_;
    vector<shared_ptr<Boundary_Source> > boundary_sources_;
};

#endif
