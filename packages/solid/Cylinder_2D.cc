#include "Cylinder_2D.hh"

#include <cmath>

#include "Vector_Functions_2D.hh"

using namespace std;

namespace vf2 = Vector_Functions_2D;

Cylinder_2D::
Cylinder_2D(int index,
            Surface_Type surface_type,
            double radius,
            vector<double> const &origin):
    Cylinder(index,
             2, // dimension
             surface_type),
    radius_(radius),
    origin_(origin),
    direction_(2, // dimension
               0)
{
}

Cylinder_2D::Relation Cylinder_2D::
relation(vector<double> const &particle_position,
         bool check_equality) const
{
    vector<double> const x = vf2::subtract(particle_position, origin_);
    
    double const r = vf2::magnitude(x);
    double const dr = r - radius_;

    if (check_equality)
    {
        if (std::abs(dr) <= relation_tolerance_)
        {
            return Relation::EQUAL;
        }
    }
    
    if (dr > 0) 
    {
        return Relation::OUTSIDE;
    }
    else // if (dr < 0)
    {
        return Relation::INSIDE;
    }
}

double Cylinder_2D::
distance(vector<double> const &position) const
{
    vector<double> const x = vf2::subtract(position, origin_);
    
    double const r = vf2::magnitude(x);
    
    return std::abs(r - radius_);
}

Cylinder_2D::Intersection Cylinder_2D::
intersection(vector<double> const &particle_position,
             vector<double> const &particle_direction) const
{
    Intersection intersection;
    
    vector<double> const k0 = vf2::subtract(particle_position,
                                            origin_);

    double const l0 = vf2::magnitude_squared(k0) - radius_ * radius_;
    double const l1 = vf2::dot(k0,
                               particle_direction);
    double const l2 = vf2::magnitude_squared(particle_direction);
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
    
    intersection.position = vf2::add(particle_position,
                                     vf2::multiply(particle_direction,
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

Cylinder_2D::Normal Cylinder_2D::
normal_direction(vector<double> const &position,
                 bool check_normal) const
{
    Normal normal;
    
    // Check if point lies on circle

    vector<double> const k0 = vf2::subtract(position,
                                            origin_);
    
    if (check_normal)
    {
        if (std::abs(vf2::magnitude_squared(k0) - radius_ * radius_) > normal_tolerance_ * radius_ * radius_)
        {
            normal.exists = false;
            return normal;
        }
    }
    
    normal.direction = vf2::normalize(k0);
    
    return normal;
}

void Cylinder_2D::
check_class_invariants() const
{
    Assert(origin_.size() == dimension_);
    Assert(direction_.size() == dimension_);
}
