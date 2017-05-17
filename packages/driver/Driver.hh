#ifndef Driver_hh
#define Driver_hh

#include <string>

/*
  Create and run transport problem from XML file
*/
class Driver
{
public:

    // Constructor
    Driver(std::string input_filename);

private:

    void run_problem();
    
    std::string input_filename_;
    std::string output_filename_;
};
    
#endif
