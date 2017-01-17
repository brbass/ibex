#include "Strong_RBF_Discretization.hh"

#include <iostream>

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "RBF_Point.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Symmetric_Sparse_Storage.hh"
#include "XML_Node.hh"

using namespace std;

Strong_RBF_Discretization::
Strong_RBF_Discretization(int dimension,
                          int number_of_points,
                          int number_of_internal_points,
                          int number_of_boundary_points,
                          int number_of_neighbors,
                          vector<int> const &internal_points,
                          vector<int> const &boundary_points,
                          vector<shared_ptr<RBF_Point> > const &rbf_points,
                          shared_ptr<Constructive_Solid_Geometry> solid_geometry):
    Spatial_Discretization(),
    dimension_(dimension),
    number_of_points_(number_of_points),
    number_of_internal_points_(number_of_internal_points),
    number_of_boundary_points_(number_of_boundary_points),
    number_of_neighbors_(number_of_neighbors),
    internal_points_(internal_points),
    boundary_points_(boundary_points),
    rbf_points_(rbf_points),
    solid_geometry_(solid_geometry),
    angular_discretization_(rbf_points[0]->material()->angular_discretization()),
    energy_discretization_(rbf_points[0]->material()->energy_discretization())
{
    check_class_invariants();
}

void Strong_RBF_Discretization::
output(XML_Node output_node) const
{
    XML_Node rbf_node = output_node.append_child("spatial_discretization");

    rbf_node.set_attribute("rbf_discretization", "type");
    rbf_node.set_child_value(dimension_, "dimension");
    rbf_node.set_child_value(number_of_points_, "number_of_points");
    rbf_node.set_child_value(number_of_boundary_points_, "number_of_boundary_points");
    rbf_node.set_child_value(number_of_internal_points_, "number_of_internal_points");
    rbf_node.set_child_vector(boundary_points_, "boundary_points", "point");
    rbf_node.set_child_vector(internal_points_, "internal_points", "point");

    XML_Node points = rbf_node.append_child("points");
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        rbf_points_[i]->output(points);
    }

    solid_geometry_->output(rbf_node);
}

void Strong_RBF_Discretization::
check_class_invariants() const
{
    Assert(number_of_points_ == number_of_boundary_points_ + number_of_internal_points_);
    Assert(internal_points_.size() == number_of_internal_points_);
    Assert(boundary_points_.size() == number_of_boundary_points_);

    for (int i = 0; i < number_of_points_; ++i)
    {
        Assert(rbf_points_[i]);
    }

    Assert(solid_geometry_);
}

