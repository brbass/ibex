#include "Plane_2D.hh"
#include "Vector_Functions_2D.hh"

#include <cmath>

using namespace std;

namespace vf2 = Vector_Functions_2D;

Plane_2D::
Plane_2D(int index,
         Surface_Type surface_type,
         vector<double> const &origin,
         vector<double> const &normal):
    Plane(index,
          2, // dimension
          surface_type),
    origin_(origin),
    normal_(normal)
{
}

Plane_2D::Relation Plane_2D::
relation(vector<double> const &particle_position,
         bool check_equality) const
{
    vector<double> const k0 = vf2::subtract(particle_position,
                                            origin_);
    double const k1 = vf2::magnitude(k0);

    double const k = vf2::dot(normal_, 
                              k0);

    if (check_equality)
    {
        if (abs(k) <= relation_tolerance_)
        {
            return Relation::EQUAL;
        }
    }
    
    if (k > 0)
    {
        return Relation::POSITIVE;
    }
    else // if (k < 0)
    {
        return Relation::NEGATIVE;
    }
}

Plane_2D::Intersection Plane_2D::
intersection(vector<double> const &particle_position,
             vector<double> const &particle_direction,
             double &distance,
             vector<double> &position) const
{
    vector<double> const k0 = vf2::subtract(origin_,
                                            particle_position);
    double const l0 = vf2::dot(k0,
                               normal_);
    double const l1 = vf2::dot(particle_direction,
                               normal_);
    
    if (abs(l1) <= intersection_tolerance_)
    {
        return Intersection::PARALLEL;
    }
    
    double const s = l0 / l1;
    
    if (s > 0)
    {
        distance = s;
    }
    else
    {
        return Intersection::NEGATIVE;
    }
    
    position = vf2::add(particle_position,
                        vf2::multiply(particle_direction,
                                      distance));
    
    return Intersection::INTERSECTS;
}

bool Plane_2D::
normal_direction(vector<double> const &position,
                 vector<double> &normal,
                 bool check_normal) const
{
    if (check_normal)
    {
        vector<double> const k0 = vf2::subtract(position, origin_);
        double const k1 = vf2::magnitude(k0);

        // Check whether point lies on line
        if (vf2::dot(normal_, k0) > normal_tolerance_ * k1)
        {
            return false;
        }
    }
    
    normal = normal_;
    
    return true;
}

void Plane_2D::
check_class_invariants() const
{
    Assert(origin_.size() == dimension_);
    Assert(normal_.size() == dimension_);
}
