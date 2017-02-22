#include "Plane_3D.hh"
#include "Vector_Functions_3D.hh"

#include <cmath>

using namespace std;

namespace vf3 = Vector_Functions_3D;

Plane_3D::
Plane_3D(int index,
         Surface_Type surface_type,
         vector<double> const &origin,
         vector<double> const &normal):
    Plane(index,
          3, // dimension
          surface_type),
    origin_(origin),
    normal_(normal)
{
}

Plane_3D::Relation Plane_3D::
relation(vector<double> const &particle_position,
         bool check_equality) const
{
    double const k = vf3::dot(normal_,
                              vf3::subtract(particle_position,
                                            origin_));

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

Plane_3D::Intersection Plane_3D::
intersection(vector<double> const &particle_position,
             vector<double> const &particle_direction) const
{
    Intersection intersection;
    
    vector<double> const k0 = vf3::subtract(origin_,
                                            particle_position);
    double const l0 = vf3::dot(k0,
                               normal_);
    double const l1 = vf3::dot(particle_direction,
                               normal_);
    
    if (abs(l1) <= intersection_tolerance_)
    {
        intersection.type = Intersection::Type::PARALLEL;
        return intersection;
    }

    double const s = l0 / l1;

    if (s > 0)
    {
        intersection.distance = s;
    }
    else
    {
        intersection.type = Intersection::Type::NEGATIVE;
        return intersection;
    }
    
    intersection.position = vf3::add(particle_position,
                                     vf3::multiply(particle_direction,
                                                   intersection.distance));
    intersection.type = Intersection::Type::INTERSECTS;
    
    return intersection;
}

Plane_3D::Normal Plane_3D::
normal_direction(vector<double> const &position,
                 bool check_normal) const
{
    Normal normal;
    
    // Check whether point lies on plane

    if (check_normal)
    {
        vector<double> const k0 = vf3::subtract(position, origin_);
        double k1 = vf3::magnitude(k0);

        if (vf3::dot(normal_, k0) > normal_tolerance_ * k1)
        {
            normal.exists = false;
            return normal;
        }
    }

    normal.direction = normal_;
    
    return normal;
}

void Plane_3D::
check_class_invariants() const
{
    Assert(origin_.size() == dimension_);
    Assert(normal_.size() == dimension_);
}
