#ifndef Material_Parser_hh
#define Material_Parser_hh

#include "Material.hh"
#include "Vector_Parser.hh"

/*
  Create a Material object from XML file
*/
class Material_Parser : public Vector_Parser<Material>
{
public:

    // Constructor
    Material_Parser(pugi::xml_node &input_file,
                    shared_ptr<Angular_Discretization> angular,
                    shared_ptr<Energy_Discretization> energy);
    
    // Return Material object
    virtual vector<shared_ptr<Material> > get_ptr() override
    {
        return materials_;
    }
    
private:
    
    vector<shared_ptr<Material> > materials_;

    shared_ptr<Angular_Discretization> angular_;
    shared_ptr<Energy_Discretization> energy_;
};

#endif
