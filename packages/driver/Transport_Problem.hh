#ifndef Transport_Problem_hh
#define Transport_Problem_hh

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "XML_Node.hh"

class Angular_Discretization;
class Energy_Discretization;
class Solid_Geometry;
class Transport_Discretization;
class Weak_RBF_Sweep;
class Weak_Spatial_Discretization;

/*
  Represents some main problem to be solved
  Lets parser pick between categories
*/
class Transport_Problem
{
public:
    
    Transport_Problem(XML_Node input_node,
                      XML_Node output_node,
                      bool print = false);

    void solve();
    
private:

    void get_weak_data(std::shared_ptr<Energy_Discretization> &energy,
                       std::shared_ptr<Angular_Discretization> &angular,
                       std::shared_ptr<Solid_Geometry> &solid,
                       std::shared_ptr<Weak_Spatial_Discretization> &spatial,
                       std::shared_ptr<Transport_Discretization> &transport,
                       std::shared_ptr<Weak_RBF_Sweep> &sweep);
    
    void solve_eigenvalue();
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
