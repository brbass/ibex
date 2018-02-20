#ifndef Manufactured_Problem_hh
#define Manufactured_Problem_hh

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "XML_Node.hh"

class Angular_Discretization;
class Energy_Discretization;
class Manufactured_Solution;
class Meshless_Sweep;
class Solid_Geometry;
class Transport_Discretization;
class Weak_Spatial_Discretization;

/*
  Represents some main problem to be solved
  Lets parser pick between categories
*/
class Manufactured_Problem
{
public:
    
    Manufactured_Problem(XML_Node input_node,
                         XML_Node output_node,
                         bool print = false);
    
    void solve();
    
private:
    
    void get_weak_data(std::shared_ptr<Energy_Discretization> &energy,
                       std::shared_ptr<Angular_Discretization> &angular,
                       std::shared_ptr<Manufactured_Solution> &solution,
                       std::shared_ptr<Solid_Geometry> &solid,
                       std::shared_ptr<Weak_Spatial_Discretization> &spatial,
                       std::shared_ptr<Transport_Discretization> &transport,
                       std::shared_ptr<Meshless_Sweep> &sweep);
    
    void solve_steady_state();
    
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
