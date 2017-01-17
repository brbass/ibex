#ifndef Driver_hh
#define Driver_hh

#include <string>
#include <vector>

/*
  Create and run transport problem from XML file
*/
class Driver
{
public:

    // Constructor
    Driver(string filename);

private:

    // Run transport problem
    void run_problem();
    
    std::string xml_in_;
    std::string xml_out_;

    void add_time(double time,
                  string description)
    {
        times_.push_back(time);
        times_description_.push_back(description);
    }
    
    std::vector<double> times_;
    std::vector<string> times_description_;
};

#endif
