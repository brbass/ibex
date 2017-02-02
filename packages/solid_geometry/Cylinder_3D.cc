#include "Cylinder_3D.hh"
#include "Vector_Functions_3D.hh"

#include <cmath>

using namespace std;

namespace vf3 = Vector_Functions_3D;

Cylinder_3D::
Cylinder_3D(int index,
            Surface_Type surface_type,
            double radius,
            vector<double> const &origin,
            vector<double> const &direction):
    Cylinder(index,
             3, // dimension
             surface_type),
    radius_(radius),
    origin_(origin),
    direction_(direction)
{
}

Cylinder_3D::Relation Cylinder_3D::
relation(vector<double> const &particle_position,
         bool check_equality) const
{
    // Cross product approach
    // vector<double> const k0 = vf3::cross(direction_,
    //                                      vf3::subtract(particle_position,
    //                                                    origin_));

    // double const r = vf3::magnitude(k0);

    // Dot product approach
    vector<double> const k0 = vf3::subtract(particle_position,
                                            origin_);
    vector<double> const n = vf3::subtract(k0,
                                           vf3::multiply(direction_,
                                                         vf3::dot(direction_,
                                                                  k0)));
    double const r = vf3::magnitude(n);
    double const dr = r - radius_;

    if (check_equality)
    {
        if (abs(dr) <= relation_tolerance_)
        {
            return Relation::EQUAL;
        }
    }
    
    if (dr > 0)
    {
        return Relation::OUTSIDE;
    }
    else // if (dr > 0)
    {
        return Relation::INSIDE;
    }
}

Cylinder_3D::Intersection Cylinder_3D::
intersection(vector<double> const &particle_position,
             vector<double> const &particle_direction,
             double &distance,
             vector<double> &position) const
{
    Intersection intersection;
    
    vector<double> const j0 = vf3::subtract(particle_position,
                                            origin_);
    vector<double> const k0 = vf3::subtract(j0,
                                            vf3::multiply(direction_,
                                                          vf3::dot(direction_,
                                                                   j0)));
    vector<double> const k1 = vf3::subtract(particle_direction,
                                            vf3::multiply(direction_,
                                                          vf3::dot(direction_,
                                                                   particle_direction)));
    
    double const l0 = vf3::magnitude_squared(k0) - radius_ * radius_;
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
                                                   distance));

    if (l3 <= intersection_tolerance_)
    {
        intersection.type = Intersection::Type::TANGEANT;
        return intersection;
    }
    else
    {
        intersection.type = Intersection::INTERSECTS;
        return intersection;
    }
}

Cylinder_3D::Normal Cylinder_3D::
normal_direction(vector<double> const &position,
                 vector<double> &normal,
                 bool check_normal) const
{
    Normal normal;
    
    vector<double> const k0 = vf3::subtract(position,
                                            origin_);
    vector<double> const n = vf3::subtract(k0,
                                           vf3::multiply(direction_,
                                                         vf3::dot(direction_,
                                                                  k0)));

    if (check_normal)
    {
        if (abs(vf3::magnitude_squared(n) - radius_ * radius_) > normal_tolerance_ * radius_ * radius_)
        {
            normal.exists = false;
            return normal;
        }
    }

    normal.exists = true;
    normal.direction = vf3::normalize(n);
    
    return normal;
}

void Cylinder_3D::
check_class_invariants() const
{
    Assert(origin_.size() == dimension_);
    Assert(direction_.size() == dimension_);
}
