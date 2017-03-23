#include "Weak_Spatial_Discretization_Factory.hh"

#include <cmath>
#include <iostream>

#include "Basis_Function.hh"
#include "Cartesian_Distance.hh"
#include "Cartesian_Plane.hh"
#include "Check.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Distance.hh"
#include "KD_Tree.hh"
#include "Linear_MLS_Function.hh"
#include "Meshless_Function_Factory.hh"
#include "RBF.hh"
#include "RBF_Factory.hh"
#include "RBF_Function.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weight_Function.hh"

using namespace std;

Weak_Spatial_Discretization_Factory::
Weak_Spatial_Discretization_Factory(shared_ptr<Constructive_Solid_Geometry> solid_geometry):
    solid_geometry_(solid_geometry)
{
    Assert(solid_geometry_->cartesian_boundaries());
}

void Weak_Spatial_Discretization_Factory::
get_basis_functions(int number_of_points,
                    vector<shared_ptr<Meshless_Function> > const &functions,
                    vector<shared_ptr<Basis_Function> > &bases) const
{
    Assert(functions.size() == number_of_points);
    
    int dimension = solid_geometry_->dimension();
    vector<shared_ptr<Cartesian_Plane> > boundaries
        = solid_geometry_->cartesian_boundary_surfaces();
    
    // Get basis functions
    bases.resize(number_of_points);
    for (int i = 0; i < number_of_points; ++i)
    {
        // Find local boundaries
        shared_ptr<Meshless_Function> function = functions[i];
        vector<shared_ptr<Cartesian_Plane> > local_boundaries;
        meshless_factory_.get_boundary_surfaces(function,
                                                boundaries,
                                                local_boundaries);
        bases[i]
            = make_shared<Basis_Function>(i,
                                          dimension,
                                          function,
                                          local_boundaries);
    }
}

void Weak_Spatial_Discretization_Factory::
get_weight_functions(int number_of_points,
                     Weight_Function::Options options,
                     vector<vector<int> > const &neighbors,
                     vector<shared_ptr<Meshless_Function> > const &functions,
                     vector<shared_ptr<Basis_Function> > const &bases,
                     vector<shared_ptr<Weight_Function> > &weights) const
{
    Assert(neighbors.size() == number_of_points);
    Assert(bases.size() == number_of_points);
    
    int dimension = solid_geometry_->dimension();
    vector<shared_ptr<Cartesian_Plane> > boundaries
        = solid_geometry_->cartesian_boundary_surfaces();

    // Get basis functions
    weights.resize(number_of_points);
    for (int i = 0; i < number_of_points; ++i)
    {
        // Find local boundaries
        shared_ptr<Meshless_Function> function = functions[i];
        vector<shared_ptr<Cartesian_Plane> > local_boundaries;
        meshless_factory_.get_boundary_surfaces(function,
                                                boundaries,
                                                local_boundaries);

        // Get local basis functions
        vector<int> const &local_neighbors = neighbors[i];
        int number_of_bases = local_neighbors.size();
        vector<shared_ptr<Basis_Function> > local_bases(number_of_bases);
        for (int j = 0; j < number_of_bases; ++j)
        {
            local_bases[j] = bases[local_neighbors[j]];
        }
        
        weights[i]
            = make_shared<Weight_Function>(i,
                                           dimension,
                                           options,
                                           function,
                                           local_bases,
                                           solid_geometry_,
                                           local_boundaries);
    }
}

shared_ptr<Weak_Spatial_Discretization> Weak_Spatial_Discretization_Factory::
get_simple_discretization(int num_dimensional_points,
                          double radius_num_intervals,
                          bool basis_mls,
                          bool weight_mls,
                          string basis_type,
                          string weight_type,
                          Weight_Function::Options options) const
{
    int dimension = solid_geometry_->dimension();
    vector<shared_ptr<Cartesian_Plane> > boundaries
        = solid_geometry_->cartesian_boundary_surfaces();
    
    // Get points
    int number_of_points;
    vector<vector<double> > points;
    vector<int> dimensional_points(dimension, num_dimensional_points);
    vector<vector<double> > limits;
    meshless_factory_.get_boundary_limits(dimension,
                                          boundaries,
                                          limits);
    meshless_factory_.get_cartesian_points(dimension,
                                           dimensional_points,
                                           limits,
                                           number_of_points,
                                           points);
    
    // Get KD tree
    shared_ptr<KD_Tree> kd_tree
        = make_shared<KD_Tree>(dimension,
                               number_of_points,
                               points);
    
    // Get neighbors
    double interval = points[1][0] - points[0][0];
    double radius = interval * radius_num_intervals;
    vector<double> radii(number_of_points, radius);
    vector<vector<int> > neighbors;
    meshless_factory_.get_neighbors(kd_tree,
                                    dimension,
                                    number_of_points,
                                    radii,
                                    radii,
                                    points,
                                    neighbors);
    
    // Get RBF and distance
    RBF_Factory rbf_factory;
    shared_ptr<RBF> basis_rbf = rbf_factory.get_rbf(basis_type); 
    shared_ptr<RBF> weight_rbf = rbf_factory.get_rbf(basis_type);
    shared_ptr<Distance> distance
        = make_shared<Cartesian_Distance>(dimension);
    
    // Get RBF functions
    vector<shared_ptr<Meshless_Function> > meshless_basis;
    if (basis_mls)
    {
        // Get simple meshless functions
        vector<shared_ptr<Meshless_Function> > simple_functions;
        meshless_factory_.get_rbf_functions(number_of_points,
                                            radii,
                                            points,
                                            basis_rbf,
                                            distance,
                                            simple_functions);
        
        // Get MLS functions
        meshless_factory_.get_mls_functions(number_of_points,
                                            simple_functions,
                                            neighbors,
                                            meshless_basis);
    }
    else
    {
        meshless_factory_.get_rbf_functions(number_of_points,
                                            radii,
                                            points,
                                            basis_rbf,
                                            distance,
                                            meshless_basis);
    }
    vector<shared_ptr<Meshless_Function> > meshless_weight;
    if (weight_mls)
    {
        // Get simple meshless functions
        vector<shared_ptr<Meshless_Function> > simple_functions;
        meshless_factory_.get_rbf_functions(number_of_points,
                                            radii,
                                            points,
                                            weight_rbf,
                                            distance,
                                            simple_functions);

        // Get MLS functions
        meshless_factory_.get_mls_functions(number_of_points,
                                            simple_functions,
                                            neighbors,
                                            meshless_weight);
    }
    else
    {
        meshless_factory_.get_rbf_functions(number_of_points,
                                            radii,
                                            points,
                                            weight_rbf,
                                            distance,
                                            meshless_weight);
    }
    
    // Get basis functions
    vector<shared_ptr<Basis_Function> > bases;
    get_basis_functions(number_of_points,
                        meshless_basis,
                        bases);

    // Get weight functions
    vector<shared_ptr<Weight_Function> > weights;
    get_weight_functions(number_of_points,
                         options,
                         neighbors,
                         meshless_weight,
                         bases,
                         weights);
    
    return make_shared<Weak_Spatial_Discretization>(bases,
                                                    weights,
                                                    kd_tree);
}
