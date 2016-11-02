#include "Point.hh"

#include "Check.hh"
#include "XML_Functions.hh"

Point::
Point(int index,
      int dimension,
      shared_ptr<Material> material,
      vector<double> const &position):
    index_(index),
    dimension_(dimension),
    point_type_(Point_Type::INTERNAL),
    material_(material),
    position_(position)
{
}

Point::
Point(int index,
      int dimension,
      shared_ptr<Material> material,
      shared_ptr<Boundary_Source> boundary_source,
      vector<double> const &position,
      vector<double> const &normal):
    index_(index),
    dimension_(dimension),
    point_type_(Point_Type::BOUNDARY),
    material_(material),
    boundary_source_(boundary_source),
    position_(position),
    normal_(normal)
{
}

void Point::
check_class_invariants() const
{
    Assert(position_.size() == dimension_);
    Assert(material_);
    
    switch(point_type_)
    {
    case Point_Type::BOUNDARY:
        Assert(normal_.size() == dimension_);
        Assert(boundary_source_);
        break;
    case Point_Type::INTERNAL:
        break;
    }
}

string Point::
point_type_string() const
{
    switch(point_type_)
    {
    case Point_Type::INTERNAL:
        return "internal";
    case Point_Type::BOUNDARY:
        return "boundary";
    }
}

void Point::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node point_node = output_node.append_child("point");

    XML_Functions::append_child(point_node, "Point", "point_type");
    XML_Functions::append_child(point_node, index_, "index");
    XML_Functions::append_child(point_node, dimension_, "dimension");
    XML_Functions::append_child(point_node, point_type_string(), "point_type");
    XML_Functions::append_child(point_node, position_, "position");
    XML_Functions::append_child(point_node, material_->index(), "material_index");
    
    switch(point_type_)
    {
    case Point_Type::BOUNDARY:
        XML_Functions::append_child(point_node, normal_, "normal");
        XML_Functions::append_child(point_node, boundary_source_->index(), "boundary_source_index");
        break;
    case Point_Type::INTERNAL:
        break;
    }
}
