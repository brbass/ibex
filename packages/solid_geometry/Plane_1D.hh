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
             vector<double> origin,
             vector<double> normal);
    
    // Plane methods
    virtual vector<double> const &normal_direction() const override
    {
        return normal_;
    }
    virtual vector<double> const &origin() const override
    {
        return origin_;
    }
    
    // Surface methods
    virtual Relation relation(vector<double> const &position,
                              bool check_equality = false) const override;
    virtual Intersection intersection(vector<double> const &initial_position,
                                      vector<double> const &initial_direction,
                                      double &distance,
                                      vector<double> &final_position) const override;
    virtual bool normal_direction(vector<double> const &position,
                                  vector<double> &normal,
                                  bool check_normal = true) const override;
    virtual void check_class_invariants() const override;
    
protected:
    
    vector<double> origin_;
    vector<double> normal_;
};

#endif
