#ifndef Angular_Discretization_Parser_hh
#define Angular_Discretization_Parser_hh

#include "Parser.hh"
#include "Angular_Discretization.hh"

/* 
   Create an object of type Angular_Discretization from xml input file
*/
class Angular_Discretization_Parser : public Parser<Angular_Discretization>
{
public:

    // Creator
    Angular_Discretization_Parser(pugi::xml_node &input_file);

    // Return pointer to created object
    virtual std::shared_ptr<Angular_Discretization> get_ptr() override
    {
        return angular_;
    }
    
private:
    
    std::shared_ptr<Angular_Discretization> angular_;
};

#endif
