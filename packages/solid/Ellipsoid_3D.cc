#include "Ellipsoid_3D.hh"

#include <cmath>

#include "Vector_Functions_3D.hh"

using namespace std;

namespace vf3 = Vector_Functions_3D;

Ellipsoid_3D::
Ellipsoid_3D(int index,
             Surface_Type surface_type,
             vector<double> const &axes,
             vector<double> const &origin):
    Ellipsoid(index,
              3, // dimension
              surface_type),
    axes_(axes),
    inverse_axes_(vf3::inverse(axes)),
    origin_(origin)
{
}

Ellipsoid_3D::Relation Ellipsoid_3D::
relation(vector<double> const &particle_position,
         bool check_equality) const
{
    vector<double> const k0
        = vf3::entrywise_product(vf3::subtract(particle_position,
                                               origin_),
                                 inverse_axes_);

    double const l0 = vf3::magnitude_squared(k0) - 1;

    if (check_equality)
    {
        if (std::abs(l0) <= relation_tolerance_)
        {
            return Relation::EQUAL;
        }
    }
    
    if (l0 > 0)
    {
        return Relation::OUTSIDE;
    }
    else // if (l0 < 0)
    {
        return Relation::INSIDE;
    }
}

Ellipsoid_3D::Intersection Ellipsoid_3D::
intersection(vector<double> const &particle_position,
             vector<double> const &particle_direction) const
{
    Intersection intersection;
    
    vector<double> const k0
        = vf3::entrywise_product(vf3::subtract(particle_position,
                                               origin_),
                                 inverse_axes_);
    vector<double> const k1
        = vf3::entrywise_product(particle_direction,
                                 inverse_axes_);

    double const l0 = vf3::magnitude_squared(k0) - 1;
    double const l1 = vf3::dot(k0,
                               k1);
    double const l2 = vf3::magnitude_squared(k1);
    double const l3 = l1 * l1 - l0 * l2;
    
    if (l3 < 0)
    {
        intersection.type = Intersection::Type::NONE;
        return intersection;
    }
    else if (l2 <= intersection_tolerance_)
    {
        intersection.type = Intersection::Type::PARALLEL;
        return intersection;
    }
    
    double const l4 = sqrt(l3);
    
    double const s1 = (-l1 + l4) / l2;
    double const s2 = (-l1 - l4) / l2;
    
    if (s2 > 0)
    {
        intersection.distance = s2;
    }
    else if (s1 > 0)
    {
        intersection.distance = s1;
    }
    else
    {
        intersection.type = Intersection::Type::NEGATIVE;
        return intersection;
    }
    
    intersection.position = vf3::add(particle_position,
                                     vf3::multiply(particle_direction,
                                                   intersection.distance));

    if (l3 <= intersection_tolerance_)
    {
        intersection.type = Intersection::Type::TANGEANT;
        return intersection;
    }
    else
    {
        intersection.type = Intersection::Type::INTERSECTS;
        return intersection;
    }
}

Ellipsoid_3D::Normal Ellipsoid_3D::
normal_direction(vector<double> const &position,
                 bool check_normal) const
{
    Normal normal;
    
    vector<double> const k0
        = vf3::entrywise_product(vf3::subtract(position,
                                               origin_),
                                 inverse_axes_);
    vector<double> const k2
        = vf3::entrywise_product(k0,
                                 inverse_axes_);
    if (check_normal)
    {
        double const l0 = vf3::magnitude_squared(k0) - 1;
        if (std::abs(l0) > normal_tolerance_)
        {
            normal.exists = false;
            return normal;
        }
    }
    
    normal.exists = true;
    normal.direction = vf3::normalize(k2);
    
    return normal;
}

void Ellipsoid_3D::
check_class_invariants() const
{
    Assert(axes_.size() == dimension_);
    Assert(inverse_axes_.size() == dimension_);
    Assert(origin_.size() == dimension_);
}
