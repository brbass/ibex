#ifndef Sphere_hh
#define Sphere_hh

#include "Surface.hh"

class Sphere : public Surface
{
public:

    Sphere(int index,
          int dimension,
          Surface_Type surface_type);
    
    // Sphere methods
    virtual double radius() const = 0;
    virtual std::vector<double> const &origin() const = 0;
    
    // Surface methods
    virtual Surface_Class surface_class() const
    {
        return Surface_Class::SPHERE;
    }
    virtual Relation relation(std::vector<double> const &position,
                              bool check_equality = false) const = 0;
    virtual Intersection intersection(std::vector<double> const &initial_position,
                                      std::vector<double> const &initial_direction,
                                      double &distance,
                                      std::vector<double> &final_position) const = 0;
    virtual bool normal_direction(std::vector<double> const &position,
                                  std::vector<double> &normal,
                                  bool check_normal = true) const = 0;
    virtual void output(XML_Node output_node) const;
};

#endif
