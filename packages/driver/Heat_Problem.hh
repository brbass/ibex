#ifndef Heat_Problem_hh
#define Heat_Problem_hh

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "XML_Node.hh"

/*
  Represents some main problem to be solved
  Lets parser pick between categories
*/
class Heat_Problem
{
public:
    
    Heat_Problem(XML_Node input_node,
                 XML_Node output_node,
                 bool print = false);

    void solve();
    
private:

    XML_Node input_node_;
    XML_Node output_node_;

    // Print
    bool print_;
    void print_message(std::string message) const;
    
    // Timing
    void output_timing();
    std::vector<std::pair<double, std::string> > times_;
};

#endif
