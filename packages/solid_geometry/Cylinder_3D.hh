#ifndef Cylinder_3D_hh
#define Cylinder_3D_hh

#include "Cylinder.hh"

/* 
   Describes an infinite cylinder with radius "radius" about the axis specified 
   by "direction" that goes through the point "origin." Two formulas are used,
   
   ||x - x0 - \Omega \cdot (x - x0) \Omega||^2 = r^2,
   ||\Omega \times (x - x0)||^2 = r^2,
   
   both of which should be valid with the magnitude applied. The first finds the
   magnitude of the vector pointing from the point x to the center of the cylinder,
   while the second finds the magnitude of the vector perpendicular to \Omega and 
   (x - x0) that runs from the origin to the surface of the cylinder. For the normal 
   direction, the first equation should be used. 
*/
class Cylinder_3D : public Cylinder
{
public:
    
    Cylinder_3D(int index,
                Surface_Type surface_type,
                double radius,
                vector<double> const &origin,
                vector<double> const &direction);

    // Cylinder methods
    virtual double radius() const override
    {
        return radius_;
    }
    virtual vector<double> const &origin() const override
    {
        return origin_;
    }
    virtual vector<double> const &direction() const override
    {
        return direction_;
    }
    
    // Surface methods
    virtual Relation relation(vector<double> const &particle_position,
                              bool check_equality = false) const override;
    virtual Intersection intersection(vector<double> const &particle_position,
                                      vector<double> const &particle_direction,
                                      double &distance,
                                      vector<double> &position) const override;
    virtual bool normal_direction(vector<double> const &position,
                                  vector<double> &normal,
                                  bool check_normal = true) const override;
    virtual void check_class_invariants() const override;
    
protected:
    
    double radius_;
    vector<double> origin_;
    vector<double> direction_;
};

#endif
