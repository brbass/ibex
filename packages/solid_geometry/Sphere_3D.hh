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
              vector<double> const &origin);

    // Sphere methods
    virtual double radius() const override
    {
        return radius_;
    }
    virtual vector<double> const &origin() const override
    {
        return origin_;
    }

    // Surface methods
    virtual Relation relation(vector<double> const &particle_position,
                              bool check_equality = false) const override;
    virtual double distance(vector<double> const &position) const override;
    virtual Intersection intersection(vector<double> const &particle_position,
                                      vector<double> const &particle_direction,
                                      double &distance,
                                      vector<double> &position) const override;
    virtual bool normal_direction(vector<double> const &position,
                                  vector<double> &normal,
                                  bool check_normal = true) const override;
    virtual void check_class_invariants() const override;

private:
    
    double radius_;
    vector<double> origin_;
};

#endif
           
