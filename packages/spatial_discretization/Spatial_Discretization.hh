#ifndef Spatial_Discretization_hh
#define Spatial_Discretization_hh

#include <memory>
#include <vector>

#include "pugixml.hh"

#include "Point.hh"

using std::shared_ptr;
using std::vector;

/*
  Pure virtual class for a spatial discretization
*/
class Spatial_Discretization
{
public:

    // Constructor
    Spatial_Discretization();
    
    // Number of spatial points to solve for in the problem
    virtual int number_of_points() const = 0;
    
    // Number of spatial points on the boundary
    virtual int number_of_boundary_points() const = 0;
    
    // Number of spatial dimensions
    virtual int dimension() const = 0;
    
    // Number of unknowns per point
    virtual int number_of_nodes() const = 0;
    
    // Spatial points on the problem boundary
    virtual vector<int> const &boundary_points() const = 0;
    
    // Return point
    virtual shared_ptr<Point> point(int point_index) const = 0;
    
    // Output data to XML file
    virtual void output(pugi::xml_node &output_node) const = 0;

    // Check class invariants
    virtual void check_class_invariants() const = 0;
};

#endif
