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
#include "Discrete_Value_Operator.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Krylov_Eigenvalue.hh"
#include "Krylov_Steady_State.hh"
#include "Linf_Convergence.hh"
#include "Material.hh"
#include "Material_Factory.hh"
#include "Material_Parser.hh"
#include "Region.hh"
#include "Solver_Factory.hh"
#include "Source_Iteration.hh"
#include "Timer.hh"
#include "Transport_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Factory.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "Weak_RBF_Sweep.hh"

namespace ce = Check_Equality;
using namespace std;

void get_pincell(bool basis_mls,
                 bool weight_mls,
                 string basis_type,
                 string weight_type,
                 shared_ptr<Weight_Function_Options> weight_options,
                 shared_ptr<Weak_Spatial_Discretization_Options> weak_options,
                 string method,
                 int dimension,
                 int angular_rule,
                 int num_dimensional_points,
                 double radius_num_intervals,
                 shared_ptr<Weak_Spatial_Discretization> &spatial,
                 shared_ptr<Angular_Discretization> &angular,
                 shared_ptr<Energy_Discretization> &energy,
                 shared_ptr<Transport_Discretization> &transport,
                 shared_ptr<Constructive_Solid_Geometry> &solid,
                 vector<shared_ptr<Material> > &materials,
                 vector<shared_ptr<Boundary_Source> > &boundary_sources,
                 shared_ptr<Weak_RBF_Sweep> &sweeper,
                 shared_ptr<Convergence_Measure> &convergence,
                 shared_ptr<Solver> &solver)
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
    materials.resize(2);
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
    // Get internal boundaries: plane for 1D, cylinder for 2D
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
    else
    {
        vector<double> origin = {0, 0};
        surfaces[2 * dimension]
            = make_shared<Cylinder_2D>(2 * dimension, // index
                                       Surface::Surface_Type::INTERNAL,
                                       length / 4, // radius
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
    else // dimension == 2
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
    solid
        = make_shared<Constructive_Solid_Geometry>(dimension,
                                                   surfaces,
                                                   regions,
                                                   materials,
                                                   boundary_sources);
    
    // Get spatial discretization
    Weak_Spatial_Discretization_Factory spatial_factory(solid,
                                                        solid->cartesian_boundary_surfaces());
    spatial = spatial_factory.get_simple_discretization(num_dimensional_points,
                                                        radius_num_intervals,
                                                        basis_mls,
                                                        weight_mls,
                                                        basis_type,
                                                        weight_type,
                                                        weight_options,
                                                        weak_options);
    
    // Get transport discretization
    transport
        = make_shared<Transport_Discretization>(spatial,
                                                angular,
                                                energy);
    
    // Get weak RBF sweep
    Weak_RBF_Sweep::Options options;
    options.solver = Weak_RBF_Sweep::Options::Solver::AMESOS;
    sweeper
        = make_shared<Weak_RBF_Sweep>(options,
                                      spatial,
                                      angular,
                                      energy,
                                      transport);

    // Get convergence method
    convergence
        = make_shared<Linf_Convergence>();
    
    // Get source iteration
    Solver_Factory solver_factory(spatial,
                                  angular,
                                  energy,
                                  transport);
    
    if (method == "krylov_eigenvalue")
    {
        solver
            = solver_factory.get_krylov_eigenvalue(sweeper);
    }
    else
    {
        AssertMsg(false, "iteration method not found");
    }
}

int test_pincell(bool basis_mls,
                 bool weight_mls,
                 string basis_type,
                 string weight_type,
                 shared_ptr<Weight_Function_Options> weight_options,
                 shared_ptr<Weak_Spatial_Discretization_Options> weak_options,
                 string method,
                 int dimension,
                 int angular_rule,
                 int num_dimensional_points,
                 double radius_num_intervals,
                 double tolerance)
{
    int checksum = 0;
    Timer timer1;
    timer1.start();
    Timer timer2;
    timer2.start();
    int w = 16;
    bool print = true;
    
    // Initialize pointers
    shared_ptr<Weak_Spatial_Discretization> spatial;
    shared_ptr<Angular_Discretization> angular;
    shared_ptr<Energy_Discretization> energy;
    shared_ptr<Transport_Discretization> transport;
    shared_ptr<Constructive_Solid_Geometry> solid;
    vector<shared_ptr<Material> > materials;
    vector<shared_ptr<Boundary_Source> > boundary_sources;
    shared_ptr<Weak_RBF_Sweep> sweeper;
    shared_ptr<Convergence_Measure> convergence;
    shared_ptr<Solver> solver;

    // Get problem
    get_pincell(basis_mls,
                weight_mls,
                basis_type,
                weight_type,
                weight_options,
                weak_options,
                method,
                dimension,
                angular_rule,
                num_dimensional_points,
                radius_num_intervals,
                spatial,
                angular,
                energy,
                transport,
                solid,
                materials,
                boundary_sources,
                sweeper,
                convergence,
                solver);
    timer2.stop();
    
    // Solve problem
    solver->solve();
    shared_ptr<Solver::Result> result
        = solver->result();
    
    // Eigenvalue problem
    double expected = dimension == 1 ? 1.13928374494 : 0.48246;
    double calculated = result->k_eigenvalue;
    if (!ce::approx(expected, calculated, tolerance))
    {
        cerr << setw(w) << "FAILED" << endl;
        checksum += 1;
    }
    timer1.stop();
    
    // Print results
    if (print)
    {
        cout << setprecision(10);
        cout << setw(w) << "init time" << setw(w) << timer2.time() << endl;
        cout << setw(w) << "total time" << setw(w) << timer1.time() << endl;
        cout << setw(w) << "expected" << setw(w) << expected << endl;
        cout << setw(w) << "eigenvalue" << setw(w) << result->k_eigenvalue << endl;
        cout << setw(w) << "pcm" << setw(w) << (expected - result->k_eigenvalue) * 1e5 << endl;
    }
    cout << endl;
    
    return checksum;
}

int run_tests(bool mls_basis,
              bool mls_weight,
              bool supg,
              string basis_type,
              string weight_type,
              double number_of_intervals)
{
    int checksum = 0;

    shared_ptr<Weight_Function_Options> weight_options
        = make_shared<Weight_Function_Options>();
    shared_ptr<Weak_Spatial_Discretization_Options> weak_options
        = make_shared<Weak_Spatial_Discretization_Options>();
    weak_options->include_supg = supg;
    
    // Test 1D eigenvalue
    {
        weak_options->integration_ordinates = 16;
        weak_options->external_integral_calculation = true;
        weight_options->tau_const = 1.0;
        weak_options->tau_scaling = Weak_Spatial_Discretization_Options::Tau_Scaling::NONE;
    
        string description = (!weak_options->include_supg
                              ? "1D, standard, "
                              : "1D, SUPG, ");
        description += basis_type + ", " + weight_type + ", ";
        description += to_string(number_of_intervals) + ", ";
        description += to_string(mls_basis) + ", ";
        description += to_string(mls_weight) + ", ";
        vector<int> point_vals = {};//4, 8, 16, 32, 64, 128, 256, 512};
        for (int number_of_points : point_vals)
        {
            cout << description << "eigenvalue, krylov, n = " << number_of_points << endl;
            checksum += test_pincell(mls_basis,
                                     mls_weight,
                                     basis_type,
                                     weight_type,
                                     weight_options,
                                     weak_options,
                                     "krylov_eigenvalue",
                                     1, // dimension
                                     256, // ordinates
                                     number_of_points, // number of points
                                     number_of_intervals, // number of intervals
                                     1e-3); // tolerance
        }
    }
    
    // Run 2D problems
    {
        weak_options->integration_ordinates = 16;
        weak_options->external_integral_calculation = true;
        weight_options->tau_const = 1.0;
        weak_options->tau_scaling = Weak_Spatial_Discretization_Options::Tau_Scaling::NONE;
        
        string description = (!weak_options->include_supg
                              ? "2D, standard, "
                              : "2D, SUPG, ");
        description += basis_type + ", " + weight_type + ", ";
        description += to_string(number_of_intervals) + ", ";
        description += to_string(mls_basis) + ", ";
        description += to_string(mls_weight) + ", ";
    
        // Test 2D eigenvalue
        vector<int> point_vals = {4, 8, 16, 32, 64, 128};
        for (int number_of_points : point_vals)
        {
            cout << description << "eigenvalue, krylov, n = " << number_of_points << endl;
            checksum += test_pincell(mls_basis,
                                     mls_weight,
                                     basis_type,
                                     weight_type,
                                     weight_options,
                                     weak_options,
                                     "krylov_eigenvalue",
                                     2, // dimension
                                     2, // angular rule
                                     number_of_points, // number of points
                                     number_of_intervals, // number of intervals
                                     1e-4); // tolerance
        }
    }
    return checksum;
}

int main(int argc, char **argv)
{
    int checksum = 0;

    if (argc != 2)
    {
        cerr << "tst_pincell [test_num]" << endl;
        return -1;
    }
    
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int test_num = atoi(argv[1]);
    switch (test_num)
    {
    case 0:
        checksum += run_tests(false,
                              false,
                              false,
                              "wendland11",
                              "wendland11",
                              3.);
        break;
    case 1:
        checksum += run_tests(false,
                              false,
                              true,
                              "wendland11",
                              "wendland11",
                              3.);
        break;
    case 2:
        checksum += run_tests(true,
                              true,
                              false,
                              "wendland11",
                              "wendland11",
                              3.);
        break;
    case 3:
        checksum += run_tests(true,
                              true,
                              true,
                              "wendland11",
                              "wendland11",
                              3.);
        break;
    }
    
    // Close MPI
    MPI_Finalize();
    
    return checksum;
}
