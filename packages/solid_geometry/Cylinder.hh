#ifndef Cylinder_hh
#define Cylinder_hh

#include "Surface.hh"

class Cylinder : public Surface
{
public:

    Cylinder(int index,
             int dimension,
             Surface_Type surface_type);
    
    // Cylinder methods
    virtual double radius() const = 0;
    virtual std::vector<double> const &origin() const = 0;
    virtual std::vector<double> const &direction() const = 0;
    
    // Surface methods
    virtual Surface_Class surface_class() const
    {
        return Surface_Class::CYLINDER;
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
    virtual void check_class_invariants() const = 0;
    virtual void output(XML_Node output_node) const;
};

#endif
