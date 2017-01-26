#include "Strong_RBF_Point.hh"

#include "Boundary_Source.hh"
#include "Material.hh"
#include "RBF_Function.hh"
#include "XML_Node.hh"

using std::shared_ptr;
using std::vector;

Strong_RBF_Point::
Strong_RBF_Point(int index,
                 int dimension,
                 int number_of_neighbors,
                 shared_ptr<Material> material,
                 shared_ptr<RBF_Function> rbf_function):
    neighbors_set_(false),
    index_(index),
    dimension_(dimension),
    number_of_neighbors_(number_of_neighbors),
    point_type_(Strong_RBF_Point::Point_Type::INTERNAL),
    material_(material),
    position_(rbf_function->position()),
    rbf_function_(rbf_function)
{
    check_class_invariants();
}

Strong_RBF_Point::
Strong_RBF_Point(int index,
                 int dimension,
                 int number_of_neighbors,
                 shared_ptr<Material> material,
                 shared_ptr<Boundary_Source> boundary_source,
                 shared_ptr<RBF_Function> rbf_function,
                 vector<double> const &normal):
    neighbors_set_(false),
    index_(index),
    dimension_(dimension),
    number_of_neighbors_(number_of_neighbors),
    point_type_(Strong_RBF_Point::Point_Type::BOUNDARY),
    material_(material),
    boundary_source_(boundary_source),
    rbf_function_(rbf_function),
    position_(rbf_function->position())
{
    check_class_invariants();
}

void Strong_RBF_Point::
set_neighbors(vector<shared_ptr<Strong_RBF_Point> > neighbor_points)
{
    neighbors_set_ = true;
    neighbor_points_ = neighbor_points;
    check_class_invariants();
}

void Strong_RBF_Point::
output(XML_Node output_node) const
{
    XML_Node point_node = output_node.append_child("point");

    point_node.set_attribute("Strong_RBF_Point", "point_type");
    point_node.set_attribute(index_, "index");
    point_node.set_child_value(dimension_, "dimension");
    point_node.set_child_value(point_type_string(), "point_type");
    point_node.set_child_vector(position_, "position");
    point_node.set_child_value(material_->index(), "material_index");

    switch(point_type_)
    {
    case Point_Type::BOUNDARY:
        point_node.set_child_vector(normal_, "normal");
        point_node.set_child_value(boundary_source_->index(), "boundary_source_index");
        break;
    case Point_Type::INTERNAL:
        break;
    }
    
    rbf_function_->output(point_node);
}

void Strong_RBF_Point::
check_class_invariants() const
{
    Assert(material_);
    Assert(rbf_function_);
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
