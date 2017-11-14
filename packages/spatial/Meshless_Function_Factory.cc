#include "Meshless_Function_Factory.hh"

#include <cmath>
#include <iostream>

#include "Cartesian_Distance.hh"
#include "Cartesian_Plane.hh"
#include "Check.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Distance.hh"
#include "KD_Tree.hh"
#include "Linear_MLS_Function.hh"
#include "Quadratic_MLS_Function.hh"
#include "RBF.hh"
#include "RBF_Factory.hh"
#include "RBF_Function.hh"

using namespace std;

Meshless_Function_Factory::
Meshless_Function_Factory()
{
}

void Meshless_Function_Factory::
get_boundary_limits(int dimension,
                    vector<shared_ptr<Cartesian_Plane> > const &boundaries,
                    vector<vector<double> > &limits) const
{
    Assert(boundaries.size() == dimension * 2);
    
    limits.assign(dimension, vector<double>(2));
    int checksum = 0;
    for (shared_ptr<Cartesian_Plane> boundary : boundaries)
    {
        int d = boundary->surface_dimension();
        double position = boundary->position();
        double normal = boundary->normal();
        int a = normal < 0 ? 0 : 1;
        limits[d][a] = position;
        
        checksum += a + 2 * d;
    }
    AssertMsg(checksum == dimension * (2 * dimension - 1),
              "boundary surface missing");
}


void Meshless_Function_Factory::
get_cartesian_points(int dimension,
                     vector<int> dimensional_points,
                     vector<vector<double> > const &limits,
                     int &number_of_points,
                     vector<vector<double> > &points) const
{
    // Check sizes
    Assert(dimensional_points.size() == dimension);
    Assert(limits.size() == dimension);
    
    // Get total number of points
    number_of_points = 1;
    for (int d = 0; d < dimension; ++d)
    {
        Assert(dimensional_points[d] > 1);
        number_of_points *= dimensional_points[d];
    }

    // Get Cartesian grid of points
    points.assign(number_of_points, vector<double>(dimension));
    int index = 0;
    vector<int> indices(dimension, 0);
    vector<double> dh(dimension);
    for (int d = 0; d < dimension; ++d)
    {
        dh[d] = (limits[d][1] - limits[d][0]) / static_cast<double>(dimensional_points[d] - 1);
    }
    while (index < number_of_points)
    {
        for (int d = 0; d < dimension; ++d)
        {
            points[index][d] = limits[d][0] + dh[d] * indices[d];
        }
        
        index += 1;
        indices[0] += 1;
        for (int d = 0; d < dimension - 1; ++d)
        {
            if (indices[d] >= dimensional_points[d])
            {
                indices[d] = 0;
                indices[d + 1] += 1;
            }
        }
    }
}

void Meshless_Function_Factory::
get_radii_nearest(shared_ptr<KD_Tree> kd_tree,
                  int dimension,
                  int number_of_points,
                  int number_of_neighbors,
                  double radius_multiplier,
                  vector<double> &radii) const
{
    radii.resize(number_of_points);
    for (int i = 0; i < number_of_points; ++i)
    {
        // Get nearest points
        vector<int> indices;
        vector<double> squared_distances;
        vector<double> const position = kd_tree->point(i);
        kd_tree->find_neighbors(number_of_neighbors,
                                position,
                                indices,
                                squared_distances);

        // Calculate radii
        radii[i] = sqrt(squared_distances[number_of_neighbors - 1]) * radius_multiplier;
    }
}

void Meshless_Function_Factory::
get_radii_coverage(shared_ptr<KD_Tree> kd_tree,
                   int dimension,
                   int number_of_points,
                   int number_of_neighbors,
                   double radius_multiplier,
                   vector<double> &radii) const
{
    radii.assign(number_of_points, 0);
    for (int i = 0; i < number_of_points; ++i)
    {
        // Get nearest points
        vector<int> indices;
        vector<double> squared_distances;
        vector<double> const position = kd_tree->point(i);
        kd_tree->find_neighbors(number_of_neighbors,
                                position,
                                indices,
                                squared_distances);
        
        for (int j = 0; j < number_of_neighbors; ++j)
        {
            int k = indices[j];
            if (radii[k] < squared_distances[j])
            {
                radii[k] = squared_distances[j];
            }
        }
    }

    for (int i = 0; i < number_of_points; ++i)
    {
        radii[i] = radius_multiplier * sqrt(radii[i]);
    }
}

void Meshless_Function_Factory::
get_neighbors(shared_ptr<KD_Tree> kd_tree,
              int dimension,
              int number_of_points,
              vector<double> const &radii,
              vector<double> const &other_radii,
              vector<vector<double> > const &positions,
              vector<vector<int> > &neighbors,
              vector<vector<double> > &squared_distances) const
{
    Assert(radii.size() == number_of_points);
    Assert(radii.size() == number_of_points);
    
    // Get maximum possible radius for neighboring points
    double max_radius = 0;
    for (double radius: other_radii)
    {
        if (radius > max_radius)
        {
            max_radius = radius;
        }
    }
    
    // Find overlapping points
    neighbors.resize(number_of_points);
    squared_distances.resize(number_of_points);
    for (int i = 0; i < number_of_points; ++i)
    {
        // Find possible overlapping points
        double const radius = radii[i];
        double search_radius = radius + max_radius;
        vector<double> const &position = positions[i];
        vector<int> local_neighbors;
        vector<double> local_squared_distances;
        kd_tree->radius_search(search_radius,
                               position,
                               local_neighbors,
                               local_squared_distances);
        Assert(local_neighbors.size() > 0);
        
        // Check to see if points actually overlap
        int number_of_local_neighbors = local_neighbors.size();
        vector<int> checked_indices;
        vector<double> checked_distances;
        checked_indices.reserve(number_of_local_neighbors);
        checked_distances.reserve(number_of_local_neighbors);
        for (int j = 0; j < number_of_local_neighbors; ++j)
        {
            int index = local_neighbors[j];
            double distance_squared = local_squared_distances[j];
            
            // Get total basis + weight radius
            double total_radius = radius + radii[index];

            // Check whether point falls inside radius
            if (total_radius * total_radius > distance_squared)
            {
                checked_indices.push_back(index);
                checked_distances.push_back(distance_squared);
            }
        }
        checked_indices.shrink_to_fit();
        checked_distances.shrink_to_fit();
        neighbors[i] = checked_indices;
        squared_distances[i] = checked_distances;
    }
}

void Meshless_Function_Factory::
get_rbf_functions(int number_of_points,
                  vector<double> const &radii,
                  vector<vector<double> > const &points,
                  shared_ptr<RBF> rbf,
                  shared_ptr<Distance> distance,
                  vector<shared_ptr<Meshless_Function> > &functions) const
{
    Assert(radii.size() == number_of_points);
    Assert(points.size() == number_of_points);

    // Get functions
    functions.resize(number_of_points);
    for (int i = 0; i < number_of_points; ++i)
    {
        // Calculate shape parameter from radii
        double const radius = radii[i];
        vector<double> const &position = points[i];
        double shape;
        switch (rbf->range())
        {
        case RBF::Range::LOCAL:
            shape = rbf->radius() / radius;
            break;
        case RBF::Range::GLOBAL:
            shape = 5. / radius;
        }
        
        // Create function
        functions[i]
            = make_shared<RBF_Function>(i,
                                        shape,
                                        position,
                                        rbf,
                                        distance);
    }
}

void Meshless_Function_Factory::
get_mls_functions(int order,
                  int number_of_points,
                  vector<shared_ptr<Meshless_Function> > const &functions,
                  vector<vector<int> > const &neighbors,
                  vector<shared_ptr<Meshless_Function> > &mls_functions) const
{
    Assert(functions.size() == number_of_points);
    Assert(neighbors.size() == number_of_points);

    // Get MLS functions
    mls_functions.resize(number_of_points);
    for (int i = 0; i < number_of_points; ++i)
    {
        // Get local neighbors
        vector<int> const &local_neighbors = neighbors[i];
        int number_of_local_neighbors = local_neighbors.size();
        vector<shared_ptr<Meshless_Function> > neighbor_functions(number_of_local_neighbors);
        for (int j = 0; j < number_of_local_neighbors; ++j)
        {
            int index = local_neighbors[j];
            neighbor_functions[j] = functions[index];
        }
        
        // Create MLS function
        switch (order)
        {
        case 1:
            mls_functions[i]
                = make_shared<Linear_MLS_Function>(neighbor_functions);
            break;
        case 2:
            mls_functions[i]
                = make_shared<Quadratic_MLS_Function>(neighbor_functions);
            break;
        default:
            AssertMsg(false, "invalid order for MLS function");
        }
    }
}

bool Meshless_Function_Factory::
check_point_conditioning(int const number_of_points,
                         vector<double> const &radii,
                         vector<vector<int> > const &neighbors,
                         vector<vector<double> > const &squared_distances) const
{
    // Check that points aren't too close together
    bool passed = true;
    double conditioning_tolerance = 1e-6;
    for (int i = 0; i < number_of_points; ++i)
    {
        double const radius = radii[i];
        double const radius2 = radius * radius;
        vector<int> const &local_neighbors = neighbors[i];
        vector<double> const &local_squared_distances = squared_distances[i];
        int const number_of_local_neighbors = local_neighbors.size();
        for (int j = 0; j < number_of_local_neighbors; ++j)
        {
            if (i != local_neighbors[j])
            {
                double const cond = local_squared_distances[j] / radius2;

                if (cond < conditioning_tolerance)
                {
                    passed = false;
                    break;
                }
            }
        }
    }

    return passed;
}

void Meshless_Function_Factory::
get_boundary_surfaces(shared_ptr<Meshless_Function> function,
                      vector<shared_ptr<Cartesian_Plane> > const &boundaries,
                      vector<shared_ptr<Cartesian_Plane> > &local_boundaries) const
{
    double const radius = function->radius();
    vector<double> const position = function->position();
    local_boundaries.clear();
    for (shared_ptr<Cartesian_Plane> boundary : boundaries)
    {
        // Get distance between boundary and center
        double distance = abs(position[boundary->surface_dimension()] - boundary->position());
        
        // If boundary is within the radius, add to local surfaces
        if (distance < radius)
        {
            local_boundaries.push_back(boundary);
        }
    }
}

void Meshless_Function_Factory::
get_cartesian_mls_functions(int order,
                            int dimension,
                            double radius_num_intervals,
                            vector<int> dimensional_points,
                            vector<vector<double> > const &limits,
                            string rbf_type,
                            vector<shared_ptr<Meshless_Function> > &functions) const
{
    // Check input data
    Assert(dimension <= 3);
    Assert(radius_num_intervals > 0);
    Assert(dimensional_points.size() == dimension);
    Assert(limits.size() == dimension);
    
    // Get points
    int number_of_points;
    vector<vector<double> > points;
    get_cartesian_points(dimension,
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
    vector<vector<double> > squared_distances;
    get_neighbors(kd_tree,
                  dimension,
                  number_of_points,
                  radii,
                  radii,
                  points,
                  neighbors,
                  squared_distances);
    
    // Get RBF and distance
    RBF_Factory rbf_factory;
    shared_ptr<RBF> rbf = rbf_factory.get_rbf(rbf_type); 
    shared_ptr<Distance> distance
        = make_shared<Cartesian_Distance>(dimension);
    
    // Get weighting functions
    vector<shared_ptr<Meshless_Function> > simple_functions;
    get_rbf_functions(number_of_points,
                      radii,
                      points,
                      rbf,
                      distance,
                      simple_functions);
    
    // Get MLS functions
    get_mls_functions(order,
                      number_of_points,
                      simple_functions,
                      neighbors,
                      functions);
}
                            
