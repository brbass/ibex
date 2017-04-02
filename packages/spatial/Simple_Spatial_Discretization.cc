#include "Simple_Spatial_Discretization.hh"

#include "Dimensional_Moments.hh"
#include "XML_Node.hh"

using std::make_shared;
using std::shared_ptr;
using std::vector;

Simple_Spatial_Discretization::
Simple_Spatial_Discretization(vector<shared_ptr<Point> > points):
    dimension_(points[0]->dimension()),
    number_of_points_(points.size()),
    points_(points)
{
    number_of_boundary_points_ = 0;
    for (int i = 0; i < number_of_points_; ++i)
    {
        switch (points_[i]->point_type())
        {
        case Point::Point_Type::INTERNAL:
            break;
        case Point::Point_Type::BOUNDARY:
            number_of_boundary_points_ += 1;
            break;
        }
    }

    dimensional_moments_
        = make_shared<Dimensional_Moments>(false, // supg
                                           dimension_);
}

void Simple_Spatial_Discretization::
output(XML_Node output_node) const
{
    XML_Node node = output_node.append_child("spatial_discretization");

    node.set_attribute("simple_spatial_discretization", "discretization_type");
    node.set_child_value(dimension_, "dimension");
    node.set_child_value(number_of_points_, "number_of_points");
    node.set_child_value(number_of_boundary_points_, "number_of_boundary_points");
    
    XML_Node points = output_node.append_child("points");
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        points_[i]->output(points);
    }
}

void Simple_Spatial_Discretization::
check_class_invariants() const
{
    Assert(points_.size() == number_of_points_);
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        Assert(points_[i]);
    }
}
