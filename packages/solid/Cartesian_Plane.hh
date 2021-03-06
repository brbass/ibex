#ifndef Cartesian_Plane_hh
#define Cartesian_Plane_hh

#include "Surface.hh"

#include <vector>

class Cartesian_Plane : public Surface
{
public:

    Cartesian_Plane(int index,
                    int dimension,
                    Surface_Type surface_type,
                    int surface_dimension,
                    double position,
                    double normal);

    // Cartesian_Plane methods
    int surface_dimension() const
    {
        return surface_dimension_;
    }
    double position() const
    {
        return position_;
    }
    double normal() const
    {
        return normal_;
    }
    
    // Surface methods
    virtual Surface_Class surface_class() const override
    {
        return Surface_Class::CARTESIAN_PLANE;
    }
    virtual Relation relation(std::vector<double> const &position,
                              bool check_equality = false) const override;
    virtual Intersection intersection(std::vector<double> const &initial_position,
                                      std::vector<double> const &initial_direction) const override;
    virtual Normal normal_direction(std::vector<double> const &position,
                                  bool check_normal = true) const override;
    virtual void check_class_invariants() const override;
    virtual void output(XML_Node output_node) const override;
    
private:
    
    int surface_dimension_;
    double position_;
    double normal_;
};

#endif
