#ifndef Plane_1D_hh
#define Plane_1D_hh

#include "Plane.hh"

/* 
   Describes a point in a 1D problem at "origin"
   with normal vector "normal".
*/
class Plane_1D : public Plane
{
public:
    
    Plane_1D(int index,
             Surface_Type surface_type,
             std::vector<double> origin,
             std::vector<double> normal);
    
    // Plane methods
    virtual std::vector<double> const &normal_direction() const override
    {
        return normal_;
    }
    virtual std::vector<double> const &origin() const override
    {
        return origin_;
    }
    
    // Surface methods
    virtual Relation relation(std::vector<double> const &position,
                              bool check_equality = false) const override;
    virtual Intersection intersection(std::vector<double> const &initial_position,
                                      std::vector<double> const &initial_direction) const override;
    virtual Normal normal_direction(std::vector<double> const &position,
                                  bool check_normal = true) const override;
    virtual void check_class_invariants() const override;
    
protected:
    
    std::vector<double> origin_;
    std::vector<double> normal_;
};

#endif
