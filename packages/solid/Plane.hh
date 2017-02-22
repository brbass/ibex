#ifndef Plane_hh
#define Plane_hh

#include <vector>

#include "Surface.hh"

class Plane : public Surface
{
public:

    Plane(int index,
          int dimension,
          Surface_Type surface_type);
    
    // Plane methods
    virtual std::vector<double> const &normal_direction() const = 0;
    virtual std::vector<double> const &origin() const = 0;
    
    // Surface methods
    virtual Surface_Class surface_class() const
    {
        return Surface_Class::PLANE;
    }
    virtual Relation relation(std::vector<double> const &position,
                              bool check_equality = false) const = 0;
    virtual Intersection intersection(std::vector<double> const &initial_position,
                                      std::vector<double> const &initial_direction) const = 0;
    virtual Normal normal_direction(std::vector<double> const &position,
                                    bool check_normal = true) const = 0;
    virtual void output(XML_Node output_node) const;
};

#endif
