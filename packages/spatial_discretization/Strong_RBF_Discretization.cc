#include "RBF_Discretization.hh"

#include <iostream>

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
// #include "Material.hh"
#include "RBF_Point.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Symmetric_Sparse_Storage.hh"
#include "XML_Functions.hh"

using std::make_shared;

RBF_Discretization::
RBF_Discretization(bool store_distances,
                   int dimension,
                   int number_of_points,
                   int number_of_internal_points,
                   int number_of_boundary_points,
                   int number_of_neighbors,
                   vector<int> const &internal_points,
                   vector<int> const &boundary_points,
                   vector<shared_ptr<RBF_Point> > const &rbf_points,
                   shared_ptr<Constructive_Solid_Geometry> solid_geometry):
    Spatial_Discretization(),
    store_distances_(store_distances),
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
    if (store_distances_)
    {
        int number_of_groups = energy_discretization_->number_of_groups();
        distances_.resize(number_of_groups);
        
        if (rbf_points[0]->rbf_function()->distance()->energy_dependent())
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                distances_[g] = make_shared<Symmetric_Sparse_Storage<int, double> >(number_of_points_);
            }
        }
        else
        {
            shared_ptr<Symmetric_Sparse_Storage<int, double> > distances = make_shared<Symmetric_Sparse_Storage<int, double> >(number_of_points_);
            for (int g = 0; g < number_of_groups; ++g)
            {
                distances_[g] = distances;
            }
        }
        // calculate_distances();
    }       
    
    check_class_invariants();
}

void RBF_Discretization::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node rbf_node = output_node.append_child("spatial_discretization");

    XML_Functions::append_child(rbf_node, "rbf_discretization", "discretization_type");
    XML_Functions::append_child(rbf_node, dimension_, "dimension");
    XML_Functions::append_child(rbf_node, number_of_points_, "number_of_points");
    XML_Functions::append_child(rbf_node, number_of_boundary_points_, "number_of_boundary_points");
    XML_Functions::append_child(rbf_node, number_of_internal_points_, "number_of_internal_points");
    XML_Functions::append_child(rbf_node, boundary_points_, "boundary_points", "point");
    XML_Functions::append_child(rbf_node, internal_points_, "internal_points", "point");

    pugi::xml_node points = rbf_node.append_child("points");
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        rbf_points_[i]->output(points);
    }

    solid_geometry_->output(rbf_node);
}

void RBF_Discretization::
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

double RBF_Discretization::
get_distance(int i,
             int i0,
             int g) const
{
    if (store_distances_)
    {
        if (!distances_[g]->element_exists(i,
                                           i0))
        {
            double dist = rbf_points_[i]->rbf_function()->distance()->distance(g,
                                                                               rbf_points_[i]->position(),
                                                                               rbf_points_[i0]->position());
            distances_[g]->add_element(i,
                                       i0,
                                       dist);
            
            return dist;
        }
        else
        {
            return distances_[g]->get_element(i, 
                                              i0);
        }
    }
    else
    {
        return rbf_points_[i]->rbf_function()->distance()->distance(g,
                                                                    rbf_points_[i]->position(),
                                                                    rbf_points_[i0]->position());
    }
}

double RBF_Discretization::
basis(int i,
      int i0,
      int g,
      int o) const
{
    shared_ptr<RBF_Point> equation = rbf_points_[i];
    shared_ptr<RBF_Point> basis = rbf_points_[i0];
    shared_ptr<RBF_Function> function = basis->rbf_function();
    vector<double> const &equation_position = equation->position();
    vector<double> const &basis_position = basis->position();
    vector<double> const &shape_parameter = basis->shape_parameter();
    vector<double> const &direction = angular_discretization_->direction(o);

    double distance = get_distance(i,
                                   i0,
                                   g);
    
    return function->basis(g,
                           shape_parameter[g],
                           equation_position,
                           basis_position,
                           direction);
}

vector<double> RBF_Discretization::
gradient_basis(int i,
               int i0,
               int g,
               int o) const
{
    shared_ptr<RBF_Point> equation = rbf_points_[i];
    shared_ptr<RBF_Point> basis = rbf_points_[i0];
    shared_ptr<RBF_Function> function = basis->rbf_function();
    vector<double> const &equation_position = equation->position();
    vector<double> const &basis_position = basis->position();
    vector<double> const &shape_parameter = basis->shape_parameter();
    vector<double> const &direction = angular_discretization_->direction(o);

    double distance = get_distance(i,
                                   i0,
                                   g);
    
    return function->gradient_basis(g,
                                    shape_parameter[g],
                                    equation_position,
                                    basis_position,
                                    direction);
}

void RBF_Discretization::
calculate_distances() const
{
    shared_ptr<RBF_Point> zeroth_point = rbf_points_[0];
    bool energy_dependent = zeroth_point->rbf_function()->distance()->energy_dependent();
    int number_of_groups = zeroth_point->material()->energy_discretization()->number_of_groups();
    shared_ptr<Distance> distance_function = zeroth_point->rbf_function()->distance();
    
    distances_.resize(number_of_groups);
    
    if (energy_dependent)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            shared_ptr<Symmetric_Sparse_Storage<int, double> > distances = make_shared<Symmetric_Sparse_Storage<int, double> >(number_of_points_);
            
            for (int i = 0; i < number_of_points_; ++i)
            {
                vector<int> neighbors = rbf_points_[i]->neighbor_indices();
                
                for (int n1 = 0; n1 < number_of_neighbors_; ++n1)
                {
                    int j = neighbors[n1];

                    vector<double> const position1 = rbf_points_[j]->position();
                    
                    for (int n2 = n1; n2 < number_of_neighbors_; ++n2)
                    {
                        int k = neighbors[n2];
                        
                        vector<double> const position2 = rbf_points_[k]->position();
                        
                        if (!distances->element_exists(j, k))
                        {
                            double value = distance_function->distance(g,
                                                                       position1,
                                                                       position2);
                            
                            distances->add_element(j, k, value);
                        }
                    }
                }
            }
            
            distances_[g] = distances;
        }
    }
    else
    {
        shared_ptr<Symmetric_Sparse_Storage<int, double> > distances = make_shared<Symmetric_Sparse_Storage<int, double> >(number_of_points_);
        
        for (int i = 0; i < number_of_points_; ++i)
        {
            vector<int> neighbors = rbf_points_[i]->neighbor_indices();
            
            for (int n1 = 0; n1 < number_of_neighbors_; ++n1)
            {
                int j = neighbors[n1];
                
                vector<double> const position1 = rbf_points_[j]->position();
                
                for (int n2 = n1; n2 < number_of_neighbors_; ++n2)
                {
                    int k = neighbors[n2];
                    
                    vector<double> const position2 = rbf_points_[k]->position();
                    
                    if (!distances->element_exists(j, k))
                    {
                        double value = distance_function->distance(0, // group
                                                                   position1,
                                                                   position2);
                        
                        distances->add_element(j, k, value);
                    }
                }
            }
        }
        
        for (int g = 0; g < number_of_groups; ++g)
        {
            distances_[g] = distances;
        }
    }
}
