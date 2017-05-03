#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Factory.hh"
#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source.hh"
#include "Boundary_Source_Parser.hh"
#include "Cartesian_Plane.hh"
#include "Check_Equality.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Constructive_Solid_Geometry_Parser.hh"
#include "Cross_Section.hh"
#include "Cylinder_2D.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Material.hh"
#include "Material_Factory.hh"
#include "Material_Parser.hh"
#include "Region.hh"
#include "Sphere_3D.hh"
#include "Timer.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Factory.hh"
#include "Weak_Spatial_Discretization_Parser.hh"

namespace ce = Check_Equality;
using namespace std;

void get_pincell(bool basis_mls,
                 bool weight_mls,
                 string basis_type,
                 string weight_type,
                 Weight_Function::Options weight_options,
                 int dimension,
                 int angular_rule,
                 int num_dimensional_points,
                 double radius_num_intervals,
                 shared_ptr<Weak_Spatial_Discretization> &spatial,
                 shared_ptr<Angular_Discretization> &angular,
                 shared_ptr<Energy_Discretization> &energy)
{
    // Set constants
    double length = 4.0;
    
    // Get angular discretization
    int number_of_moments = 1;
    Angular_Discretization_Factory angular_factory;
    angular = angular_factory.get_angular_discretization(dimension,
                                                         number_of_moments,
                                                         angular_rule);
    
    // Get energy discretization
    int number_of_groups = 1;
    energy = make_shared<Energy_Discretization>(number_of_groups);
    
    // Get material
    vector<shared_ptr<Material> > materials(2);
    Material_Factory material_factory(angular,
                                      energy);
    materials[0]
        = material_factory.get_standard_material(0, // index
                                                 {1.0}, // sigma_t
                                                 {0.84}, // sigma_s
                                                 {2.4}, // nu
                                                 {0.1}, // sigma_f
                                                 {1}, // chi
                                                 {0.0}); // internal source
    materials[1]
        = material_factory.get_standard_material(1, // index
                                                 {2.0}, // sigma_t
                                                 {1.9}, // sigma_s
                                                 {0.0}, // nu
                                                 {0.0}, // sigma_f
                                                 {0.0}, // chi
                                                 {0.0}); // internal source
    
    // Get boundary source
    Boundary_Source::Dependencies boundary_dependencies;
    vector<shared_ptr<Boundary_Source> > boundary_sources;
    if (dimension == 1)
    {
        boundary_sources.resize(2);
        boundary_sources[0]
            = make_shared<Boundary_Source>(0, // index
                                           boundary_dependencies,
                                           angular,
                                           energy,
                                           vector<double>(number_of_groups, 0), // boundary source
                                           vector<double>(number_of_groups, 1.0)); // alpha
        boundary_sources[1]
            = make_shared<Boundary_Source>(1, // index
                                           boundary_dependencies,
                                           angular,
                                           energy,
                                           vector<double>(number_of_groups, 0), // boundary source
                                           vector<double>(number_of_groups, 0.0)); // alpha
    }
    else
    {
        boundary_sources.resize(1);
        boundary_sources[0]
            = make_shared<Boundary_Source>(1, // index
                                           boundary_dependencies,
                                           angular,
                                           energy,
                                           vector<double>(number_of_groups, 0), // boundary source
                                           vector<double>(number_of_groups, 0.0)); // alpha
    }
    
    // Get solid geometry
    vector<shared_ptr<Surface> > surfaces(2 * dimension + 1);
    vector<shared_ptr<Region> > regions(2);
    // Get Cartesian boundaries
    for (int d = 0; d < dimension; ++d)
    {
        int index1 = 0 + 2 * d;
        surfaces[index1]
            = make_shared<Cartesian_Plane>(index1,
                                           dimension,
                                           Surface::Surface_Type::BOUNDARY,
                                           d,
                                           -0.5 * length,
                                           -1);
        int index2 = 1 + 2 * d;
        surfaces[index2]
            = make_shared<Cartesian_Plane>(index2,
                                           dimension,
                                           Surface::Surface_Type::BOUNDARY,
                                           d,
                                           0.5 * length,
                                           1);
    }
    
    // Get internal boundaries: plane for 1D, cylinder for 2D, sphere for 3D
    if (dimension == 1)
    {
        surfaces[2]
            = make_shared<Cartesian_Plane>(2, // index
                                           dimension,
                                           Surface::Surface_Type::INTERNAL,
                                           0, // dimension of surface
                                           0., // position of surface
                                           1.); // normal
        surfaces[0]->set_boundary_source(boundary_sources[0]);
        surfaces[1]->set_boundary_source(boundary_sources[1]);
    }
    else if (dimension == 2)
    {
        vector<double> origin = {0, 0};
        surfaces[2 * dimension]
            = make_shared<Cylinder_2D>(2 * dimension, // index
                                       Surface::Surface_Type::INTERNAL,
                                       length / 4., // radius
                                       origin);
        for (int i = 0; i < 2 * dimension; ++i)
        {
            surfaces[i]->set_boundary_source(boundary_sources[0]);
        }
    }
    else // dimension == 3
    {
        vector<double> origin = {0, 0};
        surfaces[2 * dimension]
            = make_shared<Sphere_3D>(2 * dimension, // index
                                     Surface::Surface_Type::INTERNAL,
                                     length / 4., // radius
                                     origin);
        for (int i = 0; i < 2 * dimension; ++i)
        {
            surfaces[i]->set_boundary_source(boundary_sources[0]);
        }
    }
    
    // Create regions
    if (dimension == 1)
    {
        // Fuel region
        vector<Surface::Relation> fuel_relations
            = {Surface::Relation::NEGATIVE,
               Surface::Relation::NEGATIVE};
        vector<shared_ptr<Surface> > fuel_surfaces
            = {surfaces[0],
               surfaces[2]};
        regions[0]
            = make_shared<Region>(0, // index
                                  materials[0],
                                  fuel_relations,
                                  fuel_surfaces);

        // Moderator region
        vector<Surface::Relation> mod_relations
            = {Surface::Relation::POSITIVE,
               Surface::Relation::NEGATIVE};
        vector<shared_ptr<Surface> > mod_surfaces
            = {surfaces[2],
               surfaces[1]};
        regions[1]
            = make_shared<Region>(1, // index
                                  materials[1],
                                  mod_relations,
                                  mod_surfaces);
    }
    else // (dimension == 2 || dimension == 3)
    {
        // Fuel region
        vector<Surface::Relation> fuel_relations
            = {Surface::Relation::INSIDE};
        vector<shared_ptr<Surface> > fuel_surfaces
            = {surfaces[4]};
        regions[0]
            = make_shared<Region>(0, // index
                                  materials[0],
                                  fuel_relations,
                                  fuel_surfaces);

        // Moderator region: all surfaces used
        vector<Surface::Relation> mod_relations
            = {Surface::Relation::NEGATIVE,
               Surface::Relation::NEGATIVE,
               Surface::Relation::NEGATIVE,
               Surface::Relation::NEGATIVE,
               Surface::Relation::OUTSIDE};
        regions[1]
            = make_shared<Region>(1, // index
                                  materials[1],
                                  mod_relations,
                                  surfaces);
    }
    
    // Create solid geometry
    shared_ptr<Constructive_Solid_Geometry> solid
        = make_shared<Constructive_Solid_Geometry>(dimension,
                                                   surfaces,
                                                   regions,
                                                   materials,
                                                   boundary_sources);
    
    // Get spatial discretization
    Weak_Spatial_Discretization_Factory spatial_factory(solid);
    spatial
        = spatial_factory.get_simple_discretization(num_dimensional_points,
                                                    radius_num_intervals,
                                                    basis_mls,
                                                    weight_mls,
                                                    basis_type,
                                                    weight_type,
                                                    weight_options);
}

int test_integration(bool basis_mls,
                     bool weight_mls,
                     string basis_type,
                     string weight_type,
                     Weight_Function::Options weight_options,
                     int dimension,
                     int angular_rule,
                     int num_dimensional_points,
                     double radius_num_intervals,
                     double tolerance)
{
    int checksum = 0;
    
    // Initialize pointers
    shared_ptr<Angular_Discretization> angular;
    shared_ptr<Energy_Discretization> energy;
    shared_ptr<Weak_Spatial_Discretization> spatial_internal;
    shared_ptr<Weak_Spatial_Discretization> spatial_external;

    Timer timer;
    
    // Get problem with external integration
    weight_options.external_integral_calculation = true;
    weight_options.integration_ordinates = 32;
    timer.start();
    get_pincell(basis_mls,
                weight_mls,
                basis_type,
                weight_type,
                weight_options,
                dimension,
                angular_rule,
                num_dimensional_points,
                radius_num_intervals,
                spatial_external,
                angular,
                energy);
    timer.stop();
    timer.print_time();
    
    // Get problem with internal integration
    weight_options.integration_ordinates = 32;
    weight_options.external_integral_calculation = false;
    timer.start();
    get_pincell(basis_mls,
                weight_mls,
                basis_type,
                weight_type,
                weight_options,
                dimension,
                angular_rule,
                num_dimensional_points,
                radius_num_intervals,
                spatial_internal,
                angular,
                energy);
    timer.stop();
    timer.print_time();

    shared_ptr<Weight_Function> weight_external
        = spatial_external->weight(0);
    shared_ptr<Weight_Function> weight_internal
        = spatial_internal->weight(0);

    Weight_Function::Integrals const integrals_external
        = weight_external->integrals();
    Weight_Function::Integrals const integrals_internal
        = weight_internal->integrals();

    shared_ptr<Material> material_external
        = weight_external->material();
    shared_ptr<Material> material_internal
        = weight_internal->material();
    
    int const w = 16;
    cout << "iv_b_w" << endl;
    cout << setw(w) << "external";
    cout << setw(w) << "internal";
    cout << setw(w) << "difference";
    cout << endl;
    for (int i = 0; i < integrals_external.iv_b_w.size(); ++i)
    {
        cout << setw(w) << integrals_external.iv_b_w[i];
        cout << setw(w) << integrals_internal.iv_b_w[i];
        cout << setw(w) << integrals_external.iv_b_w[i] - integrals_internal.iv_b_w[i];
        cout << endl;
    }
    
    return checksum;
}

int run_tests()
{
    int checksum = 0;

    Weight_Function::Options options;

    checksum += test_integration(true, // basis_mls
                                 true, // weight_mls
                                 "wendland11",
                                 "wendland11",
                                 options,
                                 1, // dimension
                                 2, // angular_rule
                                 9, // number_of_points
                                 3., // number_of_intervals
                                 1e-8); // tolerance
    
    return checksum;
}

int main(int argc, char **argv)
{
    int checksum = 0;
    
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Run tests
    checksum += run_tests();
    
    // Close MPI
    MPI_Finalize();
    
    return checksum;
}
