#ifndef Angular_Discretization_Factory_hh
#define Angular_Discretization_Factory_hh

#include <memory>

class Angular_Discretization;

/* 
   Create an object of type Angular_Discretization from xml input file
*/
class Angular_Discretization_Factory
{
public:
    
    // Constructor
    Angular_Discretization_Factory();
    
    // Get an angular discretization
    std::shared_ptr<Angular_Discretization> get_angular_discretization(int dimension,
                                                                       int number_of_moments,
                                                                       int angular_rule);
};

#endif
