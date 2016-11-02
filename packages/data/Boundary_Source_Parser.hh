#ifndef Boundary_Source_Parser_hh
#define Boundary_Source_Parser_hh

#include "Boundary_Source.hh"
#include "Vector_Parser.hh"

/*
  Create a Boundary_Source object from XML file
*/
class Boundary_Source_Parser : public Vector_Parser<Boundary_Source>
{
public:

    // Constructor
    Boundary_Source_Parser(pugi::xml_node &input_file,
                           shared_ptr<Angular_Discretization> angular,
                           shared_ptr<Energy_Discretization> energy);
    
    // Return pointer to object
    virtual vector<shared_ptr<Boundary_Source> > get_ptr() override
    {
        return sources_;
    }

private:
    
    vector<shared_ptr<Boundary_Source> > sources_;
    
    shared_ptr<Angular_Discretization> angular_;
    shared_ptr<Energy_Discretization> energy_;
};

#endif
