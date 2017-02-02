#include "Constructive_Solid_Geometry.hh"

#include <cmath>
#include <limits>

#include "Boundary_Source.hh"
#include "Cross_Section.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "Region.hh"
#include "String_Functions.hh"
#include "Surface.hh"
#include "Vector_Functions.hh"
#include "XML_Node.hh"

using namespace std;

namespace vf = Vector_Functions;

Constructive_Solid_Geometry::
Constructive_Solid_Geometry(int dimension,
                            vector<shared_ptr<Surface> > const &surfaces,
                            vector<shared_ptr<Region> > const &regions,
                            vector<shared_ptr<Material> > const &materials,
                            vector<shared_ptr<Boundary_Source> > const &boundary_sources):
    dimension_(dimension),
    surfaces_(surfaces),
    regions_(regions),
    materials_(materials),
    boundary_sources_(boundary_sources)
{
    optical_tolerance_ = 10 * numeric_limits<double>::epsilon();
    delta_distance_ = 1000 * numeric_limits<double>::epsilon();
    
    for (int i = 0; i < surfaces_.size(); ++i)
    {
        if (surfaces_[i]->surface_type() == Surface::Surface_Type::BOUNDARY)
        {
            boundary_surfaces_.push_back(surfaces_[i]);
        }
    }
    
    if (boundary_surfaces_.size() == 2 * dimension_)
    {
        int checksum = 0;
        for (int i = 0; i < boundary_surfaces_.size(); ++i)
        {
            if (boundary_surfaces_[i]->surface_class() != Surface::Surface_Class::PLANE)
            {
                cartesian_boundaries_ = false;
                break;
            }
            vector<double> position(dimension, 0);
            Surface::Normal normal
                = boundary_surfaces_[i]->normal_direction(position,
                                                          false);
            for (int d = 0; d < dimension_; ++d)
            {
                if (abs(abs(normal.direction[d]) - 1) < 1e-14)
                {
                    checksum += pow(2, d);
                }
            }
        }
        int checksum_actual = 2 * (pow(2, dimension_) - 1);
        if (checksum == checksum_actual)
        {
            cartesian_boundaries_ = true;
        }
        else
        {
            cartesian_boundaries_ = false;
        }
    }
    else
    {
        cartesian_boundaries_ = false;
    }
    
    check_class_invariants();
}

int Constructive_Solid_Geometry::
find_region(vector<double> const &position) const
{
    Region::Relation relation;
    
    for (int i = 0; i < regions_.size(); ++i)
    {
        relation = regions_[i]->relation(position);
        
        if (relation == Region::Relation::INSIDE)
        {
            return i;
        }
    }
    
    return NO_REGION; // no region found: outside of problem
}

int Constructive_Solid_Geometry::
find_region_including_surface(vector<double> const &position) const
{
    // Check region
    int region = find_region(position);
    if (region != NO_REGION)
    {
        return region;
    }

    // Check for surfaces
    int surface = find_surface(position);
    if (surface == NO_SURFACE)
    {
        return region;
    }
    
    // Check negative to the surface normal

    vector<double> direction;

    Surface::Normal normal
        = surfaces_[surface]->normal_direction(position,
                                               false);
    
    vector<double> check_position;
    new_position(-delta_distance_,
                 position,
                 normal.direction,
                 check_position);
    
    region = find_region(check_position);
    
    if (region != NO_REGION)
    {
        return region;
    }
    
    // Check positive to the surface normal
    new_position(delta_distance_,
                 position,
                 normal.direction,
                 check_position);
    
    region = find_region(check_position);
    
    return region;
}

int Constructive_Solid_Geometry::
find_surface(vector<double> const &position) const
{
    Surface::Relation relation;
    
    for (int i = 0; i < surfaces_.size(); ++i)
    {
        relation = surfaces_[i]->relation(position,
                                          true);
        
        if (relation == Surface::Relation::EQUAL)
        {
            return i;
        }
    }
    
    return NO_SURFACE; // particle not on a surface
}

int Constructive_Solid_Geometry::
next_intersection(int initial_region,
                  vector<double> const &initial_position,
                  vector<double> const &initial_direction,
                  int &final_region,
                  double &distance,
                  vector<double> &final_position) const
{
    if (initial_region == NO_REGION)
    {
        return next_intersection(initial_position,
                                 initial_direction,
                                 final_region,
                                 distance,
                                 final_position);
    }
    
    int best_surface = NO_SURFACE;
    distance = numeric_limits<double>::max();
    vector<double> best_position(dimension_);
    
    shared_ptr<Region> local_region = regions_[initial_region];
    int number_of_surfaces = local_region->number_of_surfaces();
    
    for (int i = 0; i < number_of_surfaces; ++i)
    {
        shared_ptr<Surface> local_surface = local_region->surface(i);

        Surface::Intersection intersection
            = local_surface->intersection(initial_position,
                                          initial_direction);
        if(intersection.type == Surface::Intersection::INTERSECTS)
        {
            if (intersection.distance < distance)
            {
                best_surface = local_surface->index();
                distance = intersection.distance;
                best_position = intersection.position;
            }
        }
    }
    
    final_position = best_position;
    
    vector<double> plus_position;
    new_position(delta_distance(),
                 final_position,
                 initial_direction,
                 plus_position);

    final_region = find_region(plus_position);
    
    return best_surface;
}

int Constructive_Solid_Geometry::
next_intersection(vector<double> const &initial_position,
                  vector<double> const &initial_direction,
                  int &final_region,
                  double &distance,
                  vector<double> &final_position) const
{
    int final_surface = NO_SURFACE;
    final_region = NO_REGION;
    distance = numeric_limits<double>::max();

    vector<double> best_position;

    for (int i = 0; i < number_of_surfaces(); ++i)
    {
        double current_distance;
        vector<double> current_position;

        Surface::Intersection intersection
            = surfaces_[i]->intersection(initial_position,
                                         initial_direction);
        if (intersection.type == Surface::Intersection::Type::INTERSECTS)
        {
            if (intersection.distance < distance)
            {
                vector<double> minus_position;
                vector<double> plus_position;
                
                new_position(-delta_distance(),
                             intersection.position,
                             initial_direction,
                             minus_position);
                new_position(delta_distance(),
                             intersection.position,
                             initial_direction,
                             plus_position);
                
                int minus_region = find_region(minus_position);
                int plus_region = find_region(plus_position);
                
                if (minus_region != plus_region)
                {
                    final_surface = i;
                    final_region = plus_region;
                    distance = intersection.distance;
                    best_position = intersection.position;
                }
            }
        }
    }

    final_position = best_position;
    
    return final_surface;
}

int Constructive_Solid_Geometry::
next_boundary(int initial_region,
              vector<double> const &initial_position,
              vector<double> const &initial_direction,
              int &boundary_region,
              double &distance,
              vector<double> &final_position) const
{
    if (initial_region == NO_REGION)
    {
        return next_intersection(initial_position,
                                 initial_direction,
                                 boundary_region,
                                 distance,
                                 final_position);
    }
    
    int surface_index = 0;
    int region_index = initial_region;
    int previous_region_index;
    vector<double> position(initial_position);
    distance = 0;

    while (surface_index != NO_SURFACE)
    {
        previous_region_index = region_index;
        
        new_position(delta_distance(),
                     position,
                     initial_direction,
                     position);
        
        distance += delta_distance();
        
        double next_distance;
        surface_index = next_intersection(region_index,
                                          position,
                                          initial_direction,
                                          region_index,
                                          next_distance,
                                          position);
        
        distance += next_distance;
        
        if (region_index == NO_REGION)
        {
            boundary_region = previous_region_index;
            final_position = position;
            
            break;
        }
    }
    
    return surface_index;
}

int Constructive_Solid_Geometry::
next_boundary(vector<double> const &initial_position,
              vector<double> const &initial_direction,
              int &boundary_region,
              double &distance,
              vector<double> &final_position) const
{
    int initial_region = find_region(initial_position);
    
    return next_boundary(initial_region,
                         initial_position,
                         initial_direction,
                         boundary_region,
                         distance,
                         final_position);
}

void Constructive_Solid_Geometry::
new_position(double distance,
             vector<double> const &initial_position,
             vector<double> const &initial_direction,
             vector<double> &final_position) const
{
    final_position.resize(dimension_);

    for (int i = 0; i < dimension_; ++i)
    {
        final_position[i] = initial_position[i] + initial_direction[i] * distance;
    }
}

void Constructive_Solid_Geometry::
optical_distance(vector<double> const &initial_position,
                 vector<double> const &final_position,
                 vector<double> &optical_distance) const
{
    // Find initial quantities
    
    vector<double> const connecting_vector = vf::subtract(final_position,
                                                          initial_position);
    double const cartesian_distance = vf::magnitude(connecting_vector);
    vector<double> const direction = vf::normalize(connecting_vector);
    
    double current_cartesian_distance = 0;
    bool point_reached = false;
    vector<double> current_position = initial_position;
    int number_of_groups = regions_[0]->material()->energy_discretization()->number_of_groups();
    
    optical_distance.assign(number_of_groups, 0);
    
    while (!point_reached)
    {
        // Keep track of the delta_distance values traversed
        
        double extra_distance = 0;
        
        // Find region
        
        new_position(delta_distance_,
                     current_position,
                     direction,
                     current_position);
        
        extra_distance += delta_distance_;
        
        int region = find_region_including_surface(current_position);
        
        if (region == NO_REGION)
        {
            string position_str;
            String_Functions::vector_to_string(position_str,
                                               current_position);
            
            AssertMsg(false, "no region found at " + position_str);
        }
        
        vector<double> const sigma_t = regions_[region]->material()->sigma_t()->data();

        // Find intersection
        
        int temp_region;
        double distance_to_surface;
        
        int next_surface = next_intersection(region,
                                             current_position,
                                             direction,
                                             temp_region,
                                             distance_to_surface,
                                             current_position);

        // Add the delta_distance values back in
        
        distance_to_surface += extra_distance;
        
        // Check if the required distance has been traversed
        // Add the distance traversed to the optical_distance
        
        if (distance_to_surface + current_cartesian_distance + optical_tolerance_ >= cartesian_distance)
        {
            double distance = cartesian_distance - current_cartesian_distance;
            
            for (int g = 0; g < number_of_groups; ++g)
            {
                optical_distance[g] += distance * sigma_t[g];
            }
            
            point_reached = true;
        }
        else
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                optical_distance[g] += distance_to_surface * sigma_t[g];
            }
            
            current_cartesian_distance += distance_to_surface;
        }
    }
}

void Constructive_Solid_Geometry::
check_class_invariants() const
{
    int number_of_surfaces = surfaces_.size();
    
    for (int i = 0; i < number_of_surfaces; ++i)
    {
        Assert(surfaces_[i]->dimension() == dimension_);
        Assert(surfaces_[i]);
    }

    int number_of_regions = regions_.size();

    for (int i = 0; i < number_of_regions; ++i)
    {
        Assert(regions_[i]);
    }
}

void Constructive_Solid_Geometry::
output(XML_Node output_node) const
{
    // Constructive solid geometry information
    
    output_node.set_child_value(dimension_, "dimension");

    // Surfaces

    XML_Node surfaces_node = output_node.append_child("surfaces");
    
    int number_of_surfaces = surfaces_.size();
    
    for (int i = 0; i < number_of_surfaces; ++i)
    {
        surfaces_[i]->output(surfaces_node.append_child("surface"));
    }

    // Regions
    
    XML_Node regions_node = output_node.append_child("regions");
    
    int number_of_regions = regions_.size();
    
    for (int i = 0; i < number_of_regions; ++i)
    {
        regions_[i]->output(regions_node.append_child("region"));
    }

    // Materials
    
    XML_Node materials_node = output_node.append_child("materials");
    
    int number_of_materials = materials_.size();
    
    for (int i = 0; i < number_of_materials; ++i)
    {
        materials_[i]->output(materials_node.append_child("material"));
    }

    // Boundary sources
    
    XML_Node boundary_sources_node = output_node.append_child("boundary_sources");
    
    int number_of_boundary_sources = boundary_sources_.size();
    
    for (int i = 0; i < number_of_boundary_sources; ++i)
    {
        boundary_sources_[i]->output(boundary_sources_node.append_child("boundary_source"));
    }
}

shared_ptr<Material> Constructive_Solid_Geometry::
material(vector<double> const &position) const
{
    return regions_[find_region(position)]->material();
}
