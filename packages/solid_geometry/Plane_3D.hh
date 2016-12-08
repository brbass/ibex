#ifndef Plane_3D_hh
#define Plane_3D_hh

#include "Plane.hh"

/* 
   Describes an infinite plane that goes through the point "origin" 
   and has surface normal "normal". The formula used is
   
   \Omega \cdot (x - x0) = 0
*/
class Plane_3D : public Plane
{
public:
    
    Plane_3D(int index,
             Surface_Type surface_type,
             std::vector<double> const &origin,
             std::vector<double> const &normal);
    
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
                                      std::vector<double> const &initial_direction,
                                      double &distance,
                                      std::vector<double> &final_position) const override;
    virtual bool normal_direction(std::vector<double> const &position,
                                  std::vector<double> &normal,
                                  bool check_normal = true) const override;
    virtual void check_class_invariants() const override;
    
protected:
    
    std::vector<double> origin_;
    std::vector<double> normal_;
};

#endif
