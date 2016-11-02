#include "Spatial_Parser_Functions.hh"

#include <algorithm>

#include "Boundary_Source.hh"
#include "Cartesian_Overlay.hh"
#include "Distance.hh"
#include "Material.hh"
#include "Optical_Distance.hh"
#include "Random_Number_Generator.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Vector_Functions_2D.hh"
#include "Vector_Functions_3D.hh"

using namespace std;
namespace vf2 = Vector_Functions_2D;
namespace vf3 = Vector_Functions_3D;

namespace // anonymous
{
    static Random_Number_Generator<double> mu_rng(-1, // lower
                                                  1, // upper
                                                  0); // seed
    static Random_Number_Generator<double> theta_rng(0, // lower
                                                     2 * M_PI, // upper
                                                     1); // seed
}

namespace Spatial_Parser_Functions
{
    void get_point(int dimension,
                   double bounding_radius,
                   vector<double> const &bounding_origin,
                   vector<double> &point)
    {
        point.resize(dimension);
        
        for (int d = 0; d < dimension; ++d)
        {
            point[d] = bounding_radius * mu_rng.scalar() + bounding_origin[d];
        }
    }
    
    void get_ray(int dimension,
                 double bounding_radius,
                 vector<double> const &bounding_origin,
                 vector<double> &origin,
                 vector<double> &direction)
        {
        double mu = mu_rng.scalar();
        double theta = theta_rng.scalar();

        vector<double> normal(dimension, 0); // outward normal of bounding sphere
        origin.assign(dimension, 0); // starting point on bounding sphere
        direction.assign(dimension, 0); // inward direction from bounding sphere
        
        switch(dimension)
        {
        case 1:
        {
            if (theta > M_PI)
            {
                origin[0] = bounding_origin[0] + bounding_radius;
                normal[0] = 1;
                direction[0] = -1;
            }
            else
            {
                origin[0] = bounding_origin[0] - bounding_radius;
                normal[0] = -1;
                direction[0] = 1;
            }
            
            break;
        }
        case 2:
        {
            origin[0] = bounding_origin[0] + bounding_radius * cos(theta);
            origin[1] = bounding_origin[1] + bounding_radius * sin(theta);
            
            normal[0] = cos(theta);
            normal[1] = sin(theta);

            direction[0] = normal[0];
            direction[1] = normal[1];

            while (vf2::dot(direction, normal) > 0)
            {
                double phi = theta_rng.scalar();
                
                direction[0] = cos(phi);
                direction[1] = sin(phi);
            }
            
            break;
        }
        case 3:
        {
            double sqrt_mu = sqrt(1 - mu * mu);
            
            origin[0] = bounding_origin[0] + bounding_radius * sqrt_mu * cos(theta);
            origin[1] = bounding_origin[1] + bounding_radius * sqrt_mu * sin(theta);
            origin[2] = bounding_origin[2] + bounding_radius * mu;

            normal[0] = sqrt_mu * cos(theta);
            normal[1] = sqrt_mu * sin(theta);
            normal[2] = mu;

            direction[0] = normal[0];
            direction[1] = normal[1];
            direction[2] = normal[2];

            while (vf3::dot(direction, normal) > 0)
            {
                double xi = mu_rng.scalar();
                double phi = theta_rng.scalar();
                double sqrt_xi = sqrt(1 - xi * xi);
                
                direction[0] = sqrt_xi * cos(phi);
                direction[1] = sqrt_xi * sin(phi);
                direction[2] = xi;
            }
            
            break;
        }
        default:
            AssertMsg(false, "dimension not found");
        }
        
    }

    void get_random_points(shared_ptr<Constructive_Solid_Geometry> solid_geometry,
                           shared_ptr<Distance> distance_metric,
                           int dimension,
                           int max_attempts,
                           double min_distance_boundary,
                           double min_distance_internal,
                           double bounding_radius,
                           vector<double> const &bounding_origin,
                           int &number_of_points,
                           int &number_of_boundary_points,
                           int &number_of_internal_points,
                           vector<int> &surfaces,
                           vector<int> &regions,
                           vector<int> &boundary_points,
                           vector<int> &internal_points,
                           vector<vector<double> > &positions,
                           vector<vector<double> > &boundary_normal,
                           vector<shared_ptr<Material> > &material,
                           vector<shared_ptr<Boundary_Source> > &boundary_source,
                           shared_ptr<Cartesian_Overlay> cartesian_overlay)
    {
        surfaces.resize(0);
        regions.resize(0);
        material.resize(0);
        boundary_points.resize(0);
        internal_points.resize(0);
        positions.resize(0);
        boundary_normal.resize(0);
        
        // Get parameters for bounding sphere
        
        vector<double> bounding_box(2 * dimension);
        for (int d = 0; d < dimension; ++d)
        {
            bounding_box[d + dimension * 0] = bounding_origin[d] - bounding_radius;
            bounding_box[d + dimension * 1] = bounding_origin[d] + bounding_radius;
        }
        
        vector<int> number_of_cells_per_dimension(dimension,
                                                  ceil(bounding_radius / min_distance_internal));
        
        cartesian_overlay = make_shared<Cartesian_Overlay>(dimension,
                                                           number_of_cells_per_dimension,
                                                           bounding_box,
                                                           distance_metric);
        
        // Find boundary points
        
        int current_boundary_point = 0;
        int current_point = 0;
        int current_internal_point = 0;
        int num_attempts = 0;
        
        while (num_attempts < max_attempts)
        {
            vector<double> position;
            vector<double> direction;
            vector<double> normal;
            
            int surface = Constructive_Solid_Geometry::NO_SURFACE;
            int boundary_region = Constructive_Solid_Geometry::NO_REGION;
            double distance;
            
            // Find ray that intersects with problem
            
            while (surface == Constructive_Solid_Geometry::NO_SURFACE)
            {
                // Get random ray in random inward-facing direction
                
                get_ray(dimension,
                        bounding_radius,
                        bounding_origin,
                        position,
                        direction);
                
                // See if ray intersects with a problem boundary
                // Move particle to the problem boundary
                
                surface = solid_geometry->next_boundary(position,
                                                        direction,
                                                        boundary_region,
                                                        distance,
                                                        position);
            }
            
            // Find all intersections of the ray
            
            while(surface != Constructive_Solid_Geometry::NO_SURFACE)
            {
                // Check whether point is too close to others

                bool has_neighbor = cartesian_overlay->has_neighbor(min_distance_boundary,
                                                                    position);
                
                if (!has_neighbor)
                {
                    surfaces.push_back(surface);
                    Check(solid_geometry->surface(surface)->normal_direction(position,
                                                                             normal));
                    material.push_back(solid_geometry->region(boundary_region)->material());
                    boundary_source.push_back(solid_geometry->surface(surface)->boundary_source());
                    boundary_points.push_back(current_point);
                    positions.push_back(position);
                    boundary_normal.push_back(normal);
                    
                    cartesian_overlay->add_point(position);
                    
                    current_boundary_point += 1;
                    current_point += 1;
                }
                else
                {
                    num_attempts += 1;
                }
                
                // Move particle ahead slightly to avoid finding the same boundary again
                
                solid_geometry->new_position(solid_geometry->delta_distance(),
                                             position,
                                             direction,
                                             position);
                
                // Find the next boundary
                
                surface = solid_geometry->next_boundary(position,
                                                        direction,
                                                        boundary_region,
                                                        distance,
                                                        position);
            }
        }
        
        // Get internal points
        
        num_attempts = 0;
        while (num_attempts < max_attempts)
        {
            int region = Constructive_Solid_Geometry::NO_REGION;
            vector<double> point;
            
            while(region == Constructive_Solid_Geometry::NO_REGION)
            {
                // Find random point

                get_point(dimension,
                          bounding_radius,
                          bounding_origin,
                          point);
                
                region = solid_geometry->find_region(point);
            }
            
            bool has_neighbor = cartesian_overlay->has_neighbor(min_distance_internal,
                                                                point);
            
            if (!has_neighbor)
            {
                regions.push_back(region);
                material.push_back(solid_geometry->region(region)->material());
                internal_points.push_back(current_point);

                positions.push_back(point);
                
                cartesian_overlay->add_point(point);
                
                current_point += 1;
                current_internal_point += 1;
            }
            else
            {
                num_attempts += 1;
            }
        }
        
        number_of_points = current_point;
        number_of_boundary_points = current_boundary_point;
        number_of_internal_points = current_internal_point;
    }
    
    
    void get_neighbor_information(int index,
                                  int dimension,
                                  int number_of_points,
                                  int number_of_groups,
                                  int number_of_neighbors,
                                  int number_to_average,
                                  shared_ptr<Distance> const distance_metric,
                                  vector<vector<double> > const &positions,
                                  vector<int> &neighbor_indices,
                                  vector<double> &mean_distance)
    {
        Assert(number_of_neighbors <= number_of_points);
        Assert(number_to_average <= number_of_neighbors);
    
        vector<int> points(number_of_points);
        vector<double> distances(number_of_points);
    
        if (distance_metric->energy_dependent())
        {
            vector<vector<double> > energy_distances(number_of_points);
        
            shared_ptr<Optical_Distance> optical_distance
                = dynamic_pointer_cast<Optical_Distance>(distance_metric);
        
            for (int i = 0; i < number_of_points; ++i)
            {
                vector<double> energy_distance = optical_distance->distance(positions[index],
                                                                            positions[i]);
            
                double sum = 0;
            
                for (int g = 0; g < number_of_groups; ++g)
                {
                    sum += energy_distance[g];
                }

                points[i] = i;
                distances[i] = sum / number_of_groups;
                energy_distances[i] = energy_distance;
            }

            partial_sort(points.begin(),
                         points.begin() + number_of_neighbors,
                         points.end(),
                         [&distances](int i1, int i2)
                         {
                             return distances[i1] < distances[i2];
                         });
    
            neighbor_indices.assign(points.begin(),
                                    points.begin() + number_of_neighbors);

            mean_distance.resize(number_of_groups);

            for (int g = 0; g < number_of_groups; ++g)
            {
                double sum = 0;
                
                for (int p = 0; p < number_to_average; ++p)
                {
                    int i = neighbor_indices[p];
                    
                    sum += energy_distances[i][g];
                }
                
                mean_distance[g] = sum / number_to_average;
            }
        }
        else
        {
            for (int i = 0; i < number_of_points; ++i)
            {
                points[i] = i;
                distances[i] = distance_metric->distance(0,
                                                         positions[index],
                                                         positions[i]);
            }

            partial_sort(points.begin(),
                         points.begin() + number_of_neighbors,
                         points.end(),
                         [&distances](int i1, int i2)
                         {
                             return distances[i1] < distances[i2];
                         });
            
            neighbor_indices.assign(points.begin(),
                                    points.begin() + number_of_neighbors);
            
            double sum = 0;
            
            for (int p = 0; p < number_to_average; ++p)
            {
                int i = neighbor_indices[p];
                
                sum += distances[i];
            }
            
            mean_distance.assign(number_of_groups, sum / number_to_average);
        }
    }

    void get_shape_parameter(int number_of_groups,
                             double shape_multiplier,
                             vector<double> &mean_distance,
                             vector<double> &shape_parameter)
    {
        shape_parameter.resize(number_of_groups);
        
        for (int g = 0; g < number_of_groups; ++g)
        {
            shape_parameter[g] = shape_multiplier / mean_distance[g];
        }
    }

    void get_optical_min_distance(int number_of_groups,
                                  double min_cartesian_distance,
                                  vector<double> const &sigma_t,
                                  vector<double> &min_distance)
    {
        min_distance.resize(number_of_groups);
        
        for (int g = 0; g < number_of_groups; ++g)
        {
            min_distance[g] = sigma_t[g] * min_cartesian_distance;
        }
    }
} // namespace Spatial_Parser_Functions
