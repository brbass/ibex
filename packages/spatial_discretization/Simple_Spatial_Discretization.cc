#include "Simple_Spatial_Discretization.hh"

#include "XML_Functions.hh"

Simple_Spatial_Discretization::
Simple_Spatial_Discretization(vector<shared_ptr<Point> > points):
    dimension_(points[0]->dimension()),
    number_of_points_(points.size()),
    points_(points)
{
    number_of_boundary_points_ = 0;
    number_of_internal_points_ = 0;
    boundary_points_.resize(0);
    internal_points_.resize(0);

    for (int i = 0; i < number_of_points_; ++i)
    {
        switch (points_[i]->point_type())
        {
        case Point::Point_Type::INTERNAL:
            number_of_internal_points_ += 1;
            internal_points_.push_back(i);
            break;
        case Point::Point_Type::BOUNDARY:
            number_of_boundary_points_ += 1;
            boundary_points_.push_back(i);
            break;
        }
    }
}

void Simple_Spatial_Discretization::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node node = output_node.append_child("spatial_discretization");

    XML_Functions::append_child(node, "simple_spatial_discretization", "discretization_type");
    XML_Functions::append_child(node, dimension_, "dimension");
    XML_Functions::append_child(node, number_of_points_, "number_of_points");
    XML_Functions::append_child(node, number_of_boundary_points_, "number_of_boundary_points");
    XML_Functions::append_child(node, number_of_internal_points_, "number_of_internal_points");
    XML_Functions::append_child(node, boundary_points_, "boundary_points", "point");
    XML_Functions::append_child(node, internal_points_, "internal_points", "point");
    
    pugi::xml_node points = output_node.append_child("points");
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        points_[i]->output(points);
    }
}

void Simple_Spatial_Discretization::
check_class_invariants() const
{
    Assert(number_of_points_ == number_of_boundary_points_ + number_of_internal_points_);
    Assert(internal_points_.size() == number_of_internal_points_);
    Assert(boundary_points_.size() == number_of_boundary_points_);
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        Assert(points_[i]);
    }
}
