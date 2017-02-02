#include "Sphere_3D.hh"

#include <cmath>

#include "Vector_Functions_3D.hh"

using namespace std;

namespace vf3 = Vector_Functions_3D;

Sphere_3D::
Sphere_3D(int index,
          Surface_Type surface_type,
          double radius,
          vector<double> const &origin):
    Sphere(index,
           3, // dimension
           surface_type),
    radius_(radius),
    origin_(origin)
{
}

Sphere_3D::Relation Sphere_3D::
relation(vector<double> const &particle_position,
         bool check_equality) const
{
    vector<double> const k0 = vf3::subtract(particle_position,
                                            origin_);
    
    double const r = vf3::magnitude(k0);
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
    else // if (dr < 0)
    {
        return Relation::INSIDE;
    }
}

double Sphere_3D::
distance(vector<double> const &position) const
{
    vector<double> const k0 = vf3::subtract(position,
                                            origin_);
    
    double const r = vf3::magnitude(k0);
    
    return abs(r - radius_);
}

Sphere_3D::Intersection Sphere_3D::
intersection(vector<double> const &particle_position,
             vector<double> const &particle_direction,
             double &distance,
             vector<double> &position) const
{
    Intersection intersection;
    
    vector<double> const k0 = vf3::subtract(particle_position,
                                            origin_);

    double const l0 = vf3::magnitude_squared(k0) - radius_ * radius_;
    double const l1 = vf3::dot(k0,
                               particle_direction);
    double const l2 = l1 * l1 - l0;

    if (l2 < 0)
    {
        intersection.type = Intersection::Type::NONE; // no intersection for line
        return intersection;
    }
    
    double const l3 = sqrt(l2);
    
    double const s1 = -l1 + l3;
    double const s2 = -l1 - l3;
    
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
        intersection.type = Intersection::Type::NEGATIVE; // intersection behind current point
    }
    
    intersection.position = vf3::add(particle_position,
                                     vf3::multiply(particle_direction,
                                                   distance));
    
    if (l2 <= intersection_tolerance_)
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

Sphere_3D::Normal Sphere_3D::
normal_direction(vector<double> const &position,
                 vector<double> &normal,
                 bool check_normal) const
{
    Normal normal;
    
    // Check if point lies on sphere
    
    vector<double> const k0 = vf3::subtract(position,
                                            origin_);

    if (check_normal)
    {
        if (abs(vf3::magnitude_squared(k0) - radius_ * radius_) > normal_tolerance_ * radius_ * radius_)
        {
            normal.exists = false;
            return normal;
        }
    }

    normal.exists = true;
    normal.direction = vf3::normalize(k0);
    
    return normal;
}

void Sphere_3D::
check_class_invariants() const
{
    Assert(origin_.size() == dimension_);
}
