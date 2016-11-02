#ifndef Cylinder_2D_hh
#define Cylinder_2D_hh

#include "Cylinder.hh"

/*
  Describes a circle with radius "radius" centered at the point "origin."
  The formula used is
  
  ||x - x0||^2 = r^2.
*/
class Cylinder_2D : public Cylinder
{
public:
    
    Cylinder_2D(int index,
                Surface_Type surface_type,
                double radius,
                vector<double> const &origin);
    
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
    virtual double distance(vector<double> const &position) const override;
    virtual void check_class_invariants() const override;
    
private:
    
    double radius_;
    vector<double> origin_;
    vector<double> direction_;
};

#endif
           
