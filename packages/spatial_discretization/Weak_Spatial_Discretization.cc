#include "Weak_Spatial_Discretization.hh"

#include "Basis_Function.hh"
#include "Check.hh"
#include "XML_Node.hh"

using namespace std;

Weak_Spatial_Discretization::
Weak_Spatial_Discretization(vector<shared_ptr<Basis_Function> > &bases,
                            vector<shared_ptr<Weight_Function> > &weights):
    number_of_points_(weights.size()),
    dimension_(weights[0]->dimension()),
    number_of_nodes_(weights[0]->number_of_nodes()),
    weights_(weights)
{
    // Get boundary weights
    int j = 0;
    boundary_weights_.resize(number_of_points_);
    for (int i = 0; i < number_of_points_; ++i)
    {
        if (weights_[i]->point_type() == Weight_Function::Point_Type::BOUNDARY)
        {
            boundary_weights_[j] = weights_[i];
            j += 1;
        }
    }
    number_of_boundary_weights_ = j;
    boundary_weights_.resize(number_of_boundary_weights_);

    // Get boundary bases
    j = 0;
    boundary_bases_.resize(number_of_points_);
    for (int i = 0; i < number_of_points_; ++i)
    {
        if (bases_[i]->number_of_boundary_surfaces() > 0)
        {
            boundary_bases_[j] = bases[i];
            j += 1;
        }
    }
    number_of_boundary_bases_ = j;
    boundary_bases_.resize(number_of_boundary_bases_);
    
    check_class_invariants();
}

void Weak_Spatial_Discretization::
check_class_invariants() const
{
    Assert(weights_.size() == number_of_points_);
    Assert(boundary_weights_.size() == number_of_boundary_weights_);
    Assert(bases_.size() == number_of_points_);
    Assert(boundary_bases_.size() == number_of_boundary_bases_);
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        Assert(weights_[i]);
        Assert(bases_[i]);
    }
    for (int i = 0; i < number_of_boundary_weights_; ++i)
    {
        Assert(boundary_weights_[i]);
    }
    for (int i = 0; i < number_of_boundary_bases_; ++i)
    {
        Assert(boundary_bases_[i]);
    }
}

void Weak_Spatial_Discretization::
output(XML_Node output_node) const
{
    output_node.set_child_value(number_of_points_, "number_of_points");
    output_node.set_child_value(number_of_boundary_weights_, "number_of_boundary_weights");
    output_node.set_child_value(number_of_boundary_bases_, "number_of_boundary_bases");
    output_node.set_child_value(dimension_, "dimension");
    output_node.set_child_value(number_of_nodes_, "number_of_nodes");
    
    XML_Node weights_node = output_node.append_child("weights");
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        weights_[i]->output(weights_node.append_child("weight"));
    }
    for (int i = 0; i < number_of_points_; ++i)
    {
        bases_[i]->output(weights_node.append_child("basis"));
    }
}
