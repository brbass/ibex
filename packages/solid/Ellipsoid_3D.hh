#ifndef Ellipsoid_3D_hh
#define Ellipsoid_3D_hh

#include "Ellipsoid.hh"

/* 
   Describes an ellipsoid centered at "origin" with semi-axes along Cartesian axes
   The formula used is
   
   ||(x - x0) \circ q||^2 = r^2,
   q = {1/a, 1/b, 1/c}
*/
class Ellipsoid_3D : public Ellipsoid
{
public:
    
    Ellipsoid_3D(int index,
                 Surface_Type surface_type,
                 std::vector<double> const &axes,
                 std::vector<double> const &origin);

    // Ellipsoid methods
    virtual std::vector<double> axes() const override
    {
        return axes_;
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
    
    std::vector<double> axes_;
    std::vector<double> inverse_axes_;
    std::vector<double> origin_;
};

#endif
           
