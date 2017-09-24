#ifndef Ellipsoid_hh
#define Ellipsoid_hh

#include "Surface.hh"

class Ellipsoid : public Surface
{
public:

    Ellipsoid(int index,
              int dimension,
              Surface_Type surface_type);
    
    // Ellipsoid methods
    virtual std::vector<double> axes() const = 0;
    virtual std::vector<double> const &origin() const = 0;
    
    // Surface methods
    virtual Surface_Class surface_class() const
    {
        return Surface_Class::ELLIPSOID;
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
