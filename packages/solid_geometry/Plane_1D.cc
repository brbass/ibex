#include "Plane_1D.hh"

#include <cmath>

using namespace std;

Plane_1D::
Plane_1D(int index,
         Surface_Type surface_type,
         vector<double> origin,
         vector<double> normal):
    Plane(index,
          1, // dimension
          surface_type),
    origin_(origin),
    normal_(normal)
{
}

Plane_1D::Relation Plane_1D::
relation(vector<double> const &particle_position,
         bool check_equality) const
{
    double const k = (particle_position[0] - origin_[0]) * normal_[0];

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

Plane_1D::Intersection Plane_1D::
intersection(vector<double> const &particle_position,
             vector<double> const &particle_direction,
             double &distance,
             vector<double> &position) const
{
    if (abs(particle_direction[0]) <= intersection_tolerance_)
    {
        return Intersection::PARALLEL;
    }
    
    distance = (origin_[0] - particle_position[0]) / particle_direction[0];
    position.assign(dimension_, origin_[0]);
    
    if (distance < 0)
    {
        return Intersection::NEGATIVE;
    }
    
    return Intersection::INTERSECTS;
}

bool Plane_1D::
normal_direction(vector<double> const &position,
                 vector<double> &normal,
                 bool check_normal) const
{
    if (check_normal)
    {
        if (abs(position[0] - origin_[0]) > normal_tolerance_)
        {
            return false;
        }
    }

    normal = normal_;
    
    return true;
}

void Plane_1D::
check_class_invariants() const
{
    Assert(origin_.size() == dimension_);
    Assert(normal_.size() == dimension_);
}
