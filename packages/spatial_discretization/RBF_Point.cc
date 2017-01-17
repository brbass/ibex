#include "RBF_Point.hh"

#include "XML_Functions.hh"

using std::shared_ptr;
using std::vector;

RBF_Point::
RBF_Point(int index,
          int dimension,
          int number_of_neighbors,
          shared_ptr<Material> material,
          shared_ptr<RBF_Function> rbf_function,
          vector<int> const &neighbor_indices,
          vector<double> const &position,
          vector<double> const &shape_parameter,
          vector<double> const &mean_distance):
    index_(index),
    dimension_(dimension),
    neighbors_set_(false),
    number_of_neighbors_(number_of_neighbors),
    material_(material),
    neighbor_indices_(neighbor_indices),
    shape_parameter_(shape_parameter),
    mean_distance_(mean_distance),
    rbf_function_(rbf_function)
{
    check_class_invariants();
}

RBF_Point::
RBF_Point(int index,
          int dimension,
          int number_of_neighbors,
          shared_ptr<Material> material,
          shared_ptr<Boundary_Source> boundary_source,
          shared_ptr<RBF_Function> rbf_function,
          vector<int> const &neighbor_indices,
          vector<double> const &position,
          vector<double> const &shape_parameter,
          vector<double> const &mean_distance,
          vector<double> const &normal):
    Point(index,
          dimension,
          material,
          boundary_source,
          position,
          normal),
    neighbors_set_(false),
    number_of_neighbors_(number_of_neighbors),
    neighbor_indices_(neighbor_indices),
    shape_parameter_(shape_parameter),
    mean_distance_(mean_distance),
    rbf_function_(rbf_function)
{
    check_class_invariants();
}

void RBF_Point::
set_neighbors(vector<shared_ptr<RBF_Point> > neighbor_points)
{
    neighbors_set_ = true;
    neighbor_points_ = neighbor_points;
}

void RBF_Point::
set_shape_multiplier(double shape_multiplier)
{
    int number_of_groups = material_->energy_discretization()->number_of_groups();
    
    for (int g = 0; g < number_of_groups; ++g)
    {
        shape_parameter_[g] = shape_multiplier / mean_distance_[g];
    }
}

void RBF_Point::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node point_node = output_node.append_child("point");

    XML_Functions::append_child(point_node, "RBF_Point", "point_type");
    XML_Functions::append_child(point_node, index_, "index");
    XML_Functions::append_child(point_node, dimension_, "dimension");
    XML_Functions::append_child(point_node, point_type_string(), "point_type");
    XML_Functions::append_child(point_node, position_, "position");
    XML_Functions::append_child(point_node, shape_parameter_, "shape_parameter");
    XML_Functions::append_child(point_node, material_->index(), "material_index");
    XML_Functions::append_child(point_node, neighbor_indices_, "neighbor_indices");

    switch(point_type_)
    {
    case Point_Type::BOUNDARY:
        XML_Functions::append_child(point_node, normal_, "normal");
        XML_Functions::append_child(point_node, boundary_source_->index(), "boundary_source_index");
        break;
    case Point_Type::INTERNAL:
        break;
    }
    
    rbf_function_->output(point_node);
}

void RBF_Point::
check_class_invariants() const
{
    Assert(material_);
    Assert(rbf_function_);
    Assert(neighbor_indices_.size() == number_of_neighbors_);
    Assert(shape_parameter_.size() == material_->energy_discretization()->number_of_groups());
    Assert(position_.size() == dimension_);
    
    switch(point_type_)
    {
    case Point_Type::BOUNDARY:
        Assert(boundary_source_);
        Assert(normal_.size() == dimension_);
        break;
    case Point_Type::INTERNAL:
        break;
    }

    if (neighbors_set_)
    {
        Assert(neighbor_points_.size() == number_of_neighbors_);

        for (int i = 0; i < number_of_neighbors_; ++i)
        {
            Assert(neighbor_points_[i]);
        }
    }

    Assert(rbf_function_);
}
