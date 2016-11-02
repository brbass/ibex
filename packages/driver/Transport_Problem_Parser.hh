#ifndef Transport_Problem_Parser_hh
#define Transport_Problem_Parser_hh

#include "Parser.hh"
#include "Transport_Problem.hh"

/*
  Parse a Transport_Problem object from XML file
*/
class Transport_Problem_Parser : public Parser<Transport_Problem>
{
public:

    // Constructor
    Transport_Problem_Parser(pugi::xml_node &input_file,
                             shared_ptr<Solver> solver);

    // Return pointer to object
    virtual shared_ptr<Transport_Problem> get_ptr() override
    {
        return transport_;
    }
    
private:
    
    shared_ptr<Transport_Problem> transport_;
};

#endif
