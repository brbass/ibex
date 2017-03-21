#include "Weak_Spatial_Discretization_Factory.hh"

#include "Basis_Function.hh"
#include "Cartesian_Distance.hh"
#include "Cartesian_Plane.hh"
#include "Check.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Distance.hh"
#include "KD_Tree.hh"
#include "Linear_MLS_Function.hh"
#include "RBF.hh"
#include "RBF_Factory.hh"
#include "RBF_Function.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weight_Function.hh"

using std::make_shared;
using std::shared_ptr;
using std::string;
using std::vector;

Weak_Spatial_Discretization_Factory::
Weak_Spatial_Discretization_Factory(shared_ptr<Constructive_Solid_Geometry> solid_geometry):
    solid_geometry_(solid_geometry)
{
    Assert(solid_geometry_->cartesian_boundaries());
}

void Weak_Spatial_Discretization_Factory::
get_cartesian_points(vector<int> dimensional_points,
                     int &number_of_points,
                     vector<vector<double> > &points) const
{
    // Get boundary surfaces
    int dimension = solid_geometry_->dimension();
    vector<shared_ptr<Cartesian_Plane> > boundaries
        = solid_geometry_->cartesian_boundary_surfaces();

    // Check sizes
    Assert(dimensional_points.size() == dimension);
    Assert(boundaries.size() == dimension * 2);

    // Get boundary limits
    vector<vector<double> > limits(dimension, vector<double>(2));
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

    // Get total number of points
    number_of_points = 1;
    for (int d = 0; d < dimension; ++d)
    {
        Assert(dimensional_points[d] > 1);
        number_of_points *= dimensional_points[d];
    }

    // Get Cartesian grid of points
    points.assign(number_of_points, vector<double>(dimension));
    switch (dimension)
    {
    case 1:
    {
        double xmin = limits[0][0];
        double xmax = limits[0][1];
        int nx = dimensional_points[0];
        for (int i = 0; i < nx; ++i)
        {
            points[i][0] = xmin + (xmax - xmin) * i / (nx - 1);
        }
        break;
    }
    case 2:
    {
        double xmin = limits[0][0];
        double xmax = limits[0][1];
        double ymin = limits[1][0];
        double ymax = limits[1][1];
        int nx = dimensional_points[0];
        int ny = dimensional_points[1];
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                int l = i + nx * j;

                points[l][0] = xmin + (xmax - xmin) * i / (nx - 1);
                points[l][1] = ymin + (ymax - ymin) * i / (ny - 1);
            }
        }
        break;
    }
    case 3:
    {
        double xmin = limits[0][0];
        double xmax = limits[0][1];
        double ymin = limits[1][0];
        double ymax = limits[1][1];
        double zmin = limits[2][0];
        double zmax = limits[2][1];
        int nx = dimensional_points[0];
        int ny = dimensional_points[1];
        int nz = dimensional_points[2];
        for (int k = 0; k < nz; ++k)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int i = 0; i < nx; ++i)
                {
                    int l = i + nx * (j + ny * k);

                    points[l][0] = xmin + (xmax - xmin) * i / (nx - 1);
                    points[l][1] = ymin + (ymax - ymin) * i / (ny - 1);
                    points[l][2] = zmin + (zmax - zmin) * i / (nz - 1);
                }
            }
        }
        break;
    }
    } // switch(dimension)
}

shared_ptr<KD_Tree> Weak_Spatial_Discretization_Factory::
get_kd_tree(int number_of_points,
            vector<vector<double> > &points) const
{
    int dimension = solid_geometry_->dimension();
    return make_shared<KD_Tree>(dimension,
                                number_of_points,
                                points);
}

void Weak_Spatial_Discretization_Factory::
get_neighbors(shared_ptr<KD_Tree> kd_tree,
              int number_of_points,
              vector<double> const &radii,
              vector<double> const &other_radii,
              vector<vector<double> > const &positions,
              vector<vector<int> > &neighbors) const
{
    Assert(radii.size() == number_of_points);
    Assert(radii.size() == number_of_points);

    int dimension = solid_geometry_->dimension();
    
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
    for (int i = 0; i < number_of_points; ++i)
    {
        // Find possible overlapping points
        double const radius = radii[i];
        double search_radius = radius + max_radius;
        vector<double> const &position = positions[i];
        vector<int> indices;
        vector<double> distances;
        kd_tree->radius_search(search_radius,
                               position,
                               indices,
                               distances);
        Assert(indices.size() > 0);
        
        // Check to see if points actually overlap
        vector<int> checked_indices;
        checked_indices.reserve(indices.size());
        for (int index : indices)
        {
            // Get distance squared between centers
            double distance_squared = 0;
            for (int d = 0; d < dimension; ++d)
            {
                double delta = positions[i][d] - positions[index][d];
                distance_squared += delta * delta;
            }

            // Get total basis + weight radius
            double total_radius = radius + radii[index];

            // Check whether point falls inside radius
            if (total_radius * total_radius > distance_squared)
            {
                checked_indices.push_back(index);
            }
        }
        checked_indices.shrink_to_fit();
        neighbors[i] = checked_indices;
    }
}

void Weak_Spatial_Discretization_Factory::
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
        double const shape = rbf->radius() / radius;

        // Create function
        functions[i]
            = make_shared<RBF_Function>(shape,
                                        position,
                                        rbf,
                                        distance);
    }
}

void Weak_Spatial_Discretization_Factory::
get_mls_functions(int number_of_points,
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
        vector<shared_ptr<Meshless_Function> > neighbor_functions(local_neighbors.size());
        for (int index : local_neighbors)
        {
            neighbor_functions[i] = functions[index];
        }

        // Create MLS function
        mls_functions[i]
            = make_shared<Linear_MLS_Function>(neighbor_functions);
    }
}

void Weak_Spatial_Discretization_Factory::
get_boundary_surfaces(shared_ptr<Meshless_Function> function,
                      vector<shared_ptr<Cartesian_Plane> > &local_boundaries) const
{
    vector<shared_ptr<Cartesian_Plane> > boundaries
        = solid_geometry_->cartesian_boundary_surfaces();
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

void Weak_Spatial_Discretization_Factory::
get_basis_functions(int number_of_points,
                    vector<shared_ptr<Meshless_Function> > const &functions,
                    vector<shared_ptr<Basis_Function> > &bases) const
{
    Assert(functions.size() == number_of_points);
    
    int dimension = solid_geometry_->dimension();

    // Get basis functions
    bases.resize(number_of_points);
    for (int i = 0; i < number_of_points; ++i)
    {
        // Find local boundaries
        shared_ptr<Meshless_Function> function = functions[i];
        vector<shared_ptr<Cartesian_Plane> > local_boundaries;
        get_boundary_surfaces(function,
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

    // Get basis functions
    weights.resize(number_of_points);
    for (int i = 0; i < number_of_points; ++i)
    {
        // Find local boundaries
        shared_ptr<Meshless_Function> function = functions[i];
        vector<shared_ptr<Cartesian_Plane> > local_boundaries;
        get_boundary_surfaces(function,
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
    
    // Get points
    int number_of_points;
    vector<vector<double> > points;
    vector<int> dimensional_points(dimension, num_dimensional_points);
    get_cartesian_points(dimensional_points,
                         number_of_points,
                         points);
    
    // Get KD tree
    shared_ptr<KD_Tree> kd_tree = get_kd_tree(number_of_points,
                                              points);
    
    // Get neighbors
    double interval = points[1][0] - points[0][0];
    double radius = interval * radius_num_intervals;
    vector<double> radii(number_of_points, radius);
    vector<vector<int> > neighbors;
    get_neighbors(kd_tree,
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
        get_rbf_functions(number_of_points,
                          radii,
                          points,
                          basis_rbf,
                          distance,
                          simple_functions);

        // Get MLS functions
        get_mls_functions(number_of_points,
                          simple_functions,
                          neighbors,
                          meshless_basis);
    }
    else
    {
        get_rbf_functions(number_of_points,
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
        get_rbf_functions(number_of_points,
                          radii,
                          points,
                          weight_rbf,
                          distance,
                          simple_functions);

        // Get MLS functions
        get_mls_functions(number_of_points,
                          simple_functions,
                          neighbors,
                          meshless_weight);
    }
    else
    {
        get_rbf_functions(number_of_points,
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
