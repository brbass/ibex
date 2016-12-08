#include "Plane.hh"

#include "Boundary_Source.hh"
#include "XML_Node.hh"

Plane::
Plane(int index,
      int dimension,
      Surface_Type surface_type):
    Surface(index,
            dimension,
            surface_type)
{
}

void Plane::
output(XML_Node output_node) const
{
    output_node.set_attribute(index_, "index");
    output_node.set_attribute("plane", "shape");
    output_node.set_child_value(dimension_, "dimension");
    output_node.set_child_vector(normal_direction(), "normal_direction");
    output_node.set_child_vector(origin(), "origin");

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
