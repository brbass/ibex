#ifndef Sphere_3D_hh
#define Sphere_3D_hh

#include "Sphere.hh"

/* 
   Describes a sphere centered at "origin" with radius "radius".
   The formula used is
   
   ||x - x0||^2 = r^2. 
*/
class Sphere_3D : public Sphere
{
public:
    
    Sphere_3D(int index,
              Surface_Type surface_type,
              double radius,
              std::vector<double> const &origin);

    // Sphere methods
    virtual double radius() const override
    {
        return radius_;
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

private:
    
    double radius_;
    std::vector<double> origin_;
};

#endif
           
