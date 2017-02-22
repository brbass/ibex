#include "Check_Equality.hh"
#include "Cylinder_2D.hh"
#include "Cylinder_3D.hh"
#include "Plane_1D.hh"
#include "Plane_2D.hh"
#include "Plane_3D.hh"
#include "Sphere_3D.hh"
#include "Surface.hh"
#include "Vector_Functions.hh"

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>

namespace ce = Check_Equality;
namespace vf = Vector_Functions;

using namespace std;

double tolerance()
{
    return 100 * numeric_limits<double>::epsilon();
}
double delta_distance()
{
    return 100 * numeric_limits<double>::epsilon();
}

// Make shapes such that {1, 1, 1} is on the surface,
// {1 + delta_distance(), 1, 1} is outside the surface and
// {1 - delta_distance(), 1, 1} is inside the surface.
// The normal at {1, 1, 1} should be {1, 0, 0}

shared_ptr<Surface> get_cylinder_2d()
{
    int index = 0;
    Surface::Surface_Type surface_type = Surface::Surface_Type::INTERNAL;
    double radius = 2.;
    vector<double> origin = {-1., 1.};
    
    return make_shared<Cylinder_2D>(index,
                                    surface_type,
                                    radius,
                                    origin);

}

shared_ptr<Surface> get_cylinder_3d()
{
    int index = 0;
    Surface::Surface_Type surface_type = Surface::Surface_Type::INTERNAL;
    double radius = 2.;
    vector<double> origin = {-1., 1., 7.};
    vector<double> direction = {0., 0., 1.};
    
    return make_shared<Cylinder_3D>(index,
                                    surface_type,
                                    radius,
                                    origin,
                                    direction);
}

shared_ptr<Surface> get_plane_1d()
{
    int index = 0;
    Surface::Surface_Type surface_type = Surface::Surface_Type::INTERNAL;
    vector<double> origin = {1};
    vector<double> normal = {1};

    return make_shared<Plane_1D>(index,
                                 surface_type,
                                 origin,
                                 normal);
}

shared_ptr<Surface> get_plane_2d()
{
    int index = 0;
    Surface::Surface_Type surface_type = Surface::Surface_Type::INTERNAL;
    vector<double> origin = {1., 22.};
    vector<double> normal = {1., 0.};

    return make_shared<Plane_2D>(index,
                                 surface_type,
                                 origin,
                                 normal);
}

shared_ptr<Surface> get_plane_3d()
{
    int index = 0;
    Surface::Surface_Type surface_type = Surface::Surface_Type::INTERNAL;
    vector<double> origin = {1., 22., -77.};
    vector<double> normal = {1., 0., 0.};

    return make_shared<Plane_3D>(index,
                                 surface_type,
                                 origin,
                                 normal);
}

shared_ptr<Surface> get_sphere_3d()
{
    int index = 0;
    Surface::Surface_Type surface_type = Surface::Surface_Type::INTERNAL;
    double radius = 4.;
    vector<double> origin = {-3., 1., 1.};
    
    return make_shared<Sphere_3D>(index,
                                  surface_type,
                                  radius,
                                  origin);
}

int test_surface(shared_ptr<Surface> surface,
                 string description)
{
    int checksum = 0;
    
    int dimension = surface->dimension();
    vector<double> equal_position(dimension, 1);
    vector<double> normal_direction(dimension, 0);
    normal_direction[0] = 1.;
    
    // Test relation

    {
        vector<double> position = equal_position;
        if (surface->relation(position,
                              true)
            != Surface::Relation::EQUAL)
        {
            cout << description << " EQUAL failed" << endl;
            checksum += 1;
        }

        position = equal_position;
        position[0] += delta_distance();
        if (surface->relation(position,
                              false)
            != Surface::Relation::OUTSIDE)
        {
            cout << description << " OUTSIDE failed" << endl;
            checksum += 1;
        }

        position = equal_position;
        position[0] -= delta_distance();
        if (surface->relation(position,
                              false)
            != Surface::Relation::INSIDE)
        {
            cout << description << " INSIDE failed" << endl;
            checksum += 1;
        }
    }
    
    // Test intersection and reflection
    
    {
        // Intersection

        vector<double> initial_position(dimension, 0);
        initial_position[0] = 2;
        vector<double> initial_direction
            = vf::normalize(vf::subtract(equal_position,
                                         initial_position));
        
        Surface::Intersection const intersection
            = surface->intersection(initial_position,
                                    initial_direction);
        if (intersection.type != Surface::Intersection::Type::INTERSECTS)
        {
            cout << description << " intersection failed" << endl;
            checksum += 1;
        }
        if (!ce::approx(intersection.distance, sqrt(dimension), tolerance()))
        {
            cout << description << " distance incorrect: " << intersection.distance << endl;
            checksum += 1;
        }

        if (!ce::approx(intersection.position, equal_position, tolerance()))
        {
            cout << description << " intersection position incorrect" << endl;
            cout << intersection.position[0] << endl;
            checksum += 1;
        }

        // Reflection
        
        Surface::Reflection const reflection
            = surface->reflected_direction(equal_position,
                                           initial_direction,
                                           true) ;
        if (!reflection.exists)
        {
            cout << description << " reflected normal not found" << endl;
            checksum += 1;
        }

        vector<double> expected_final_direction = initial_direction;
        expected_final_direction[0] = -expected_final_direction[0];

        if (!ce::approx(expected_final_direction, reflection.direction, tolerance()))
        {
            cout << description << " reflected direction incorrect" << endl;
            checksum += 1;
        }
    }

    // Test normal direction
    
    {
        Surface::Normal const normal
            = surface->normal_direction(equal_position,
                                        true);

        if (!normal.exists)
        {
            cout << description << " normal surface not found" << endl;
            checksum += 1;
        }

        if (!ce::approx(normal.direction, normal_direction, tolerance()))
        {
            cout << description << " normal direction incorrect" << endl;
            checksum += 1;
        }
    }
                                      
    return checksum;
}

int main()
{
    int checksum = 0;
    
    checksum += test_surface(get_cylinder_2d(), "cylinder_2d");
    checksum += test_surface(get_cylinder_3d(), "cylinder_3d");
    checksum += test_surface(get_plane_1d(), "plane_1d");
    checksum += test_surface(get_plane_2d(), "plane_2d");
    checksum += test_surface(get_plane_3d(), "plane_3d");
    checksum += test_surface(get_sphere_3d(), "plane_3d");

    return checksum;
}
