#include "Cartesian_Plane.hh"

#include <cmath>
#include <vector>

#include "Boundary_Source.hh"
#include "XML_Node.hh"

using namespace std;

Cartesian_Plane::
Cartesian_Plane(int index,
                int dimension,
                Surface_Type surface_type,
                int surface_dimension,
                double position,
                double normal):
    Surface(index,
            dimension,
            surface_type),
    surface_dimension_(surface_dimension),
    position_(position),
    normal_(normal)
{
}

Cartesian_Plane::Relation Cartesian_Plane::
relation(vector<double> const &position,
         bool check_equality) const
{
    double const k = (position[surface_dimension_] - position_) * normal_;

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
    else
    {
        return Relation::NEGATIVE;
    }
}

Cartesian_Plane::Intersection Cartesian_Plane::
intersection(vector<double> const &initial_position,
             vector<double> const &initial_direction,
             double &distance,
             vector<double> &final_position) const
{
    if (abs(initial_direction[surface_dimension_]) <= intersection_tolerance_)
    {
        return Intersection::PARALLEL;
    }

    distance = (position_ - initial_position[surface_dimension_]) / initial_direction[surface_dimension_];

    final_position.assign(dimension_, 0);
    final_position[surface_dimension_] = position_;

    if (distance < 0)
    {
        return Intersection::NEGATIVE;
    }

    return Intersection::INTERSECTS;
}

bool Cartesian_Plane::
normal_direction(vector<double> const &position,
                 vector<double> &normal,
                 bool check_normal) const
{
    if (check_normal)
    {
        if (abs(position[surface_dimension_] - position_) > normal_tolerance_)
        {
            return false;
        }
    }
    
    normal.assign(dimension_, 0);
    normal[surface_dimension_] = 1.;

    return true;
}

void Cartesian_Plane::
output(XML_Node output_node) const
{
    output_node.set_attribute(index_, "index");
    output_node.set_attribute("cartesian_plane", "shape");
    output_node.set_child_value(dimension_, "dimension");
    output_node.set_child_value(surface_dimension_, "surface_dimension");
    output_node.set_child_value(position_, "position");
    output_node.set_child_value(normal_, "normal");

    switch(surface_type_)
    {
    case Surface_Type::BOUNDARY:
        output_node.set_attribute("boundary", "type");
        output_node.set_child_value(boundary_source_->index(), "boundary_source_index");
        break;
    case Surface_Type::INTERNAL:
        output_node.set_attribute("internal", "type");
        break;
    }
}

void Cartesian_Plane::
check_class_invariants() const
{
    Assert(surface_dimension_ < dimension_);
    Assert(boundary_source_);
}
