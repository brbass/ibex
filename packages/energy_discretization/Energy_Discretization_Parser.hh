#ifndef Energy_Discretization_Parser_hh
#define Energy_Discretization_Parser_hh

#include "Energy_Discretization.hh"
#include "Parser.hh"

/*
  Create an Energy_Discretization object from XML input file
*/
class Energy_Discretization_Parser : public Parser<Energy_Discretization>
{
public:

    // Constructor
    Energy_Discretization_Parser(pugi::xml_node &input_file);

    // Return pointer to Energy_Discretization object
    virtual std::shared_ptr<Energy_Discretization> get_ptr() override
    {
        return energy_;
    }
    
private:
    
    std::shared_ptr<Energy_Discretization> energy_;
};

#endif
