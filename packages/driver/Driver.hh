#ifndef Driver_hh
#define Driver_hh

#include <string>
#include <vector>

using std::string;
using std::vector;

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
    
    string xml_in_;
    string xml_out_;

    void add_time(double time,
                  string description)
    {
        times_.push_back(time);
        times_description_.push_back(description);
    }
    
    vector<double> times_;
    vector<string> times_description_;
};

#endif
