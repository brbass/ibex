#include "Basis_Function.hh"

#include "Cartesian_Plane.hh"
#include "Conversion.hh"
#include "Meshless_Function.hh"
#include "XML_Node.hh"

using namespace std;

Basis_Function::
Basis_Function(int index,
               int dimension,
               shared_ptr<Meshless_Function> meshless_function,
               vector<shared_ptr<Cartesian_Plane> > boundary_surfaces):
    index_(index),
    boundary_index_(-1),
    dimension_(dimension),
    number_of_boundary_surfaces_(boundary_surfaces.size()),
    radius_(meshless_function->radius()),
    position_(meshless_function->position()),
    meshless_function_(meshless_function),
    boundary_surfaces_(boundary_surfaces)
{
    point_type_ = (number_of_boundary_surfaces_ == 0
                   ? Point_Type::INTERNAL
                   : Point_Type::BOUNDARY);
    check_class_invariants();
}

void Basis_Function::
check_class_invariants() const
{
    Assert(number_of_boundary_surfaces_ == boundary_surfaces_.size());
    Assert(position_.size() == dimension_);
    Assert(meshless_function_);
    for (int i = 0; i < number_of_boundary_surfaces_; ++i)
    {
        Assert(boundary_surfaces_[i]);
    }
}

void Basis_Function::
output(XML_Node output_node) const
{
    output_node.set_attribute(index_, "index");
    output_node.set_child_value(dimension_, "dimension");
    output_node.set_child_value(number_of_boundary_surfaces_, "number_of_boundary_surfaces");
    output_node.set_child_value(radius_, "radius");
    output_node.set_child_vector(position_, "position");
    meshless_function_->output(output_node.append_child("function"));

    vector<int> boundary_surfaces(number_of_boundary_surfaces_);
    for (int i = 0; i < number_of_boundary_surfaces_; ++i)
    {
        boundary_surfaces[i] = boundary_surfaces_[i]->index();
    }
    output_node.set_child_vector(boundary_surfaces, "boundary_surfaces");
}

shared_ptr<Conversion<Basis_Function::Point_Type, string> > Basis_Function::
point_type_conversion() const
{
    vector<pair<Point_Type, string> > conversions
        = {{Point_Type::BOUNDARY, "boundary"},
           {Point_Type::INTERNAL, "internal"}};
    return make_shared<Conversion<Point_Type, string> >(conversions);
}
