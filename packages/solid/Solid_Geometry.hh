#ifndef Solid_Geometry_hh
#define Solid_Geometry_hh

#include <memory>
#include <vector>

class Material;
class XML_Node;

class Solid_Geometry
{
public:

    enum Geometry_Errors
    {
        NO_REGION = -1,
        NO_SURFACE = -2
    };
    
    Solid_Geometry();

    virtual int dimension() const = 0;
    virtual int find_region(std::vector<double> const &position) const = 0;
    virtual int find_surface(std::vector<double> const &position) const = 0;
    virtual int next_intersection(std::vector<double> const &initial_position,
                                  std::vector<double> const &initial_direction,
                                  int &final_region,
                                  double &distance,
                                  std::vector<double> &final_position) const = 0;
    virtual int next_boundary(std::vector<double> const &initial_position,
                              std::vector<double> const &initial_direction,
                              int &boundary_region,
                              double &distance,
                              std::vector<double> &final_position) const = 0;
    virtual void optical_distance(std::vector<double> const &initial_position,
                                  std::vector<double> const &final_position,
                                  std::vector<double> &optical_distance) const = 0;
    virtual std::shared_ptr<Material> material(std::vector<double> const &position) const = 0;
    virtual void check_class_invariants() const = 0;
    virtual void output(XML_Node output_node) const = 0;
};

#endif
