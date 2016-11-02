#include "Sphere.hh"

Sphere::
Sphere(int index,
       int dimension,
       Surface_Type surface_type):
    Surface(index,
            dimension,
            surface_type)
{
}

void Sphere::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node node = output_node.append_child("surface");

    XML_Functions::append_child(node, index_, "index");
    XML_Functions::append_child(node, dimension_, "dimension");
    XML_Functions::append_child(node, "sphere", "type");
    XML_Functions::append_child(node, radius(), "radius");
    XML_Functions::append_child(node, origin(), "origin");

    switch(surface_type_)
    {
    case Surface_Type::BOUNDARY:
        XML_Functions::append_child(node, "boundary", "surface_type");
        XML_Functions::append_child(node, boundary_source_->index(), "boundary_source_index");
        break;
    case Surface_Type::INTERNAL:
        XML_Functions::append_child(node, "internal", "surface_type");
        break;
    }
}
