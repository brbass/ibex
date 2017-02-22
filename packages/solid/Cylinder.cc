#include "Cylinder.hh"

#include "Boundary_Source.hh"
#include "XML_Node.hh"

Cylinder::
Cylinder(int index,
         int dimension,
         Surface_Type surface_type):
    Surface(index,
            dimension,
            surface_type)
{
}

void Cylinder::
output(XML_Node output_node) const
{
    output_node.set_attribute(index_, "index");
    output_node.set_attribute("cylinder", "shape");
    output_node.set_child_value(dimension_, "dimension");
    output_node.set_child_value(radius(), "radius");
    output_node.set_child_vector(origin(), "origin");
    output_node.set_child_vector(direction(), "direction");

    switch(surface_type_)
    {
    case Surface_Type::BOUNDARY:
        output_node.set_attribute("boundary", "type");
        output_node.set_child_value(boundary_source_->index(), "boundary_source");
        break;
    case Surface_Type::INTERNAL:
        output_node.set_attribute("internal", "type");
        break;
    }
}
