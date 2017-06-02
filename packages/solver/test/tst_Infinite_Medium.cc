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
#include "Transport_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Factory.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "Weak_RBF_Sweep.hh"

namespace ce = Check_Equality;
using namespace std;

void get_one_region(bool basis_mls,
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
                    double sigma_t,
                    double sigma_s,
                    double chi_nu_sigma_f,
                    double internal_source,
                    double boundary_source,
                    double alpha,
                    double length,
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
    materials.resize(1);
    Material_Factory material_factory(angular,
                                      energy);
    materials[0]
        = material_factory.get_standard_material(0, // index
                                                 {sigma_t},
                                                 {sigma_s}, // sigma_s
                                                 {1}, // nu
                                                 {chi_nu_sigma_f}, // sigma_f
                                                 {1}, // chi
                                                 {internal_source});
    
    // Get boundary source
    boundary_sources.resize(1);
    Boundary_Source::Dependencies boundary_dependencies;
    boundary_sources[0]
        = make_shared<Boundary_Source>(0, // index
                                       boundary_dependencies,
                                       angular,
                                       energy,
                                       vector<double>(number_of_groups, boundary_source), // boundary source
                                       vector<double>(number_of_groups, alpha)); // alpha
    
    // Get solid geometry
   vector<shared_ptr<Surface> > surfaces(2 * dimension);
    vector<shared_ptr<Region> > regions(1);
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
    for (shared_ptr<Surface> surface : surfaces)
    {
        surface->set_boundary_source(boundary_sources[0]);
    }
    vector<Surface::Relation> surface_relations(2 * dimension,
                                                Surface::Relation::NEGATIVE);
    regions[0]
        = make_shared<Region>(0, // index
                              materials[0],
                              surface_relations,
                              surfaces);
    solid
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
                                                    weight_options,
                                                    weak_options);
    
    // Get transport discretization
    transport
        = make_shared<Transport_Discretization>(spatial,
                                                angular,
                                                energy);
    
    // Get weak RBF sweep
    Weak_RBF_Sweep::Options sweep_options;
    sweeper
        = make_shared<Weak_RBF_Sweep>(sweep_options,
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
    
    if (method == "krylov_steady_state")
    {
        solver
            = solver_factory.get_krylov_steady_state(sweeper,
                                                     convergence);
    }
    else if (method == "source_iteration")
    {
        solver
            = solver_factory.get_source_iteration(sweeper,
                                                  convergence);
    }
    else if (method == "krylov_eigenvalue")
    {
        solver
            = solver_factory.get_krylov_eigenvalue(sweeper);
    }
    else
    {
        AssertMsg(false, "iteration method not found");
    }
}

int test_infinite(bool basis_mls,
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
                  double sigma_t,
                  double sigma_s,
                  double chi_nu_sigma_f,
                  double internal_source,
                  double boundary_source,
                  double alpha,
                  double length,
                  double tolerance)
{
    int checksum = 0;
    
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
    get_one_region(basis_mls,
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
                   sigma_t,
                   sigma_s,
                   chi_nu_sigma_f,
                   internal_source,
                   boundary_source,
                   alpha,
                   length,
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

    // Solve problem
    solver->solve();
    shared_ptr<Solver::Result> result
        = solver->result();
    
    // Print and check results
    bool print = true;
    int w = 16;
    
    // Eigenvalue problem
    if (method.find("eigenvalue") != string::npos)
    {
        // Check that eigenvalue is correct
        double expected = chi_nu_sigma_f / (sigma_t - sigma_s);
        double calculated = result->k_eigenvalue;
        if (!ce::approx(expected, calculated, tolerance))
        {
            cerr << "eigenvalue incorrect" << endl;
            checksum += 1;
        }

        // Print results
        if (print)
        {
            cout << setw(w) << "expected" << setw(w) << expected << endl;
            cout << setw(w) << "eigenvalue" << setw(w) << result->k_eigenvalue << endl;
        }
        cout << endl;
    }
    // Steady state problem
    else
    {
        int phi_size = transport->phi_size();
        int num_values = result->phi.size();

        // Check values
        double solution = internal_source / (sigma_t - sigma_s - chi_nu_sigma_f);
        vector<double> solution_vec(phi_size, solution);

        for (int j = 0; j < num_values; ++j)
        {
            if (!ce::approx(solution_vec, result->phi[j], tolerance))
            {
                checksum += 1;
                cerr << "flux incorrect" << endl;
            }
        }
        
        // Print results
        if (print)
        {
            cout << setw(w) << "expected: " << setw(w) << solution << endl;
            cout << setw(w) << "value" << setw(w) << "weighted" << endl;
            for (int i = 0; i < phi_size; ++i)
            {
                for (int j = 0; j < num_values; ++j)
                {
                    cout << setw(w) << result->phi[j][i];
                }
                cout << endl;
            }
        }
        cout << endl;
    }
    
    return checksum;
}

int run_tests()
{
    int checksum = 0;
    
    // Run 1D problems for regular and SUPG options
    for (int i = 0; i < 2; ++i)
    {
        shared_ptr<Weight_Function_Options> weight_options
            = make_shared<Weight_Function_Options>();
        shared_ptr<Weak_Spatial_Discretization_Options> weak_options
            = make_shared<Weak_Spatial_Discretization_Options>();
        weak_options->tau_scaling
            = (i == 0
               ? Weak_Spatial_Discretization_Options::Output::STANDARD
               : Weak_Spatial_Discretization_Options::Output::SUPG);
        weak_options->integration_ordinates = 16;
        weight_options->tau_const = 1.0;
        weak_options->tau_scaling = Weak_Spatial_Discretization_Options::Tau_Scaling::NONE;
        
        string description = (weak_options->output == Weak_Spatial_Discretization_Options::Output::STANDARD
                              ? "1D standard "
                              : "1D SUPG ");
        
        // Test 1D eigenvalue
        cout << description << "eigenvalue, krylov" << endl;
        checksum += test_infinite(true, // mls basis
                                  true, // mls weight
                                  "wendland11", // basis type
                                  "wendland11", // weight type
                                  weight_options,
                                  weak_options
                                  "krylov_eigenvalue",
                                  1, // dimension
                                  16, // ordinates
                                  5, // number of points
                                  3, // number of intervals
                                  2.0, // sigma_t
                                  0.8, // sigma_s
                                  1.1, // nu_sigma_f
                                  0.0, // internal source
                                  0.0, // boundary source
                                  1.0, // alpha
                                  2.0, // length
                                  1e-4); // tolerance

        // Test 1D steady state with reflecting boundaries
        cout << description << "steady state with reflecting boundaries, krylov" << endl;
        checksum += test_infinite(true, // mls basis
                                  true, // mls weight
                                  "wendland11", // basis type
                                  "wendland11", // weight type
                                  weight_options,
                                  weak_options,
                                  "krylov_steady_state",
                                  1, // dimension
                                  16, // ordinates
                                  5, // number of points
                                  3, // number of intervals
                                  2.0, // sigma_t
                                  0.8, // sigma_s
                                  1.1, // nu_sigma_f
                                  1.0, // internal source
                                  0.0, // boundary source
                                  1.0, // alpha
                                  2.0, // length
                                  1e-4); // tolerance

        // Test 1D steady state with boundary source
        cout << description << "steady state with boundary source, source iteration" << endl;
        checksum += test_infinite(true, // mls basis
                                  true, // mls weight
                                  "wendland11", // basis type
                                  "wendland11", // weight type
                                  weight_options,
                                  weak_options,
                                  "source_iteration",
                                  1, // dimension
                                  16, // ordinates
                                  5, // number of points
                                  3, // number of intervals
                                  2.0, // sigma_t
                                  0.8, // sigma_s
                                  1.1, // nu_sigma_f
                                  1.0, // internal source
                                  1.0 / (2 * (2.0 - 0.8 - 1.1)), // boundary source
                                  0.0, // alpha
                                  2.0, // length
                                  1e-4); // tolerance
    }

    // Run 2D problems
    for (int i = 0; i < 2; ++i)
    {
        shared_ptr<Weight_Function_Options> weight_options
            = make_shared<Weight_Function_Options>();
        shared_ptr<Weak_Spatial_Discretization_Options> weak_options
            = make_shared<Weak_Spatial_Discretization_Options>();
        weak_options->tau_scaling
            = (i == 0
               ? Weak_Spatial_Discretization_Options::Output::STANDARD
               : Weak_Spatial_Discretization_Options::Output::SUPG);
        weak_options->integration_ordinates = 32;
        weight_options->tau_const = 1.0;
        weak_options->tau_scaling = Weak_Spatial_Discretization_Options::Tau_Scaling::NONE;
        
        string description = (weak_options->output == Weak_Spatial_Discretization_Options::Output::STANDARD
                              ? "2D standard "
                              : "2D SUPG ");
        
        // Test 2D eigenvalue
        cout << description << "eigenvalue, krylov" << endl;
        checksum += test_infinite(true, // mls basis
                                  true, // mls weight
                                  "wendland11", // basis type
                                  "wendland11", // weight type
                                  weight_options,
                                  "krylov_eigenvalue",
                                  2, // dimension
                                  3, // angular rule
                                  5, // number of points
                                  3, // number of intervals
                                  2.0, // sigma_t
                                  0.8, // sigma_s
                                  1.1, // nu_sigma_f
                                  0.0, // internal source
                                  0.0, // boundary source
                                  1.0, // alpha
                                  2.0, // length
                                  1e-4); // tolerance

        // Test 2D steady state with reflecting boundaries
        cout << description << "steady state with reflecting boundaries, krylov" << endl;
        checksum += test_infinite(true, // mls basis
                                  true, // mls weight
                                  "wendland11", // basis type
                                  "wendland11", // weight type
                                  weight_options,
                                  "krylov_steady_state",
                                  2, // dimension
                                  3, // angular rule
                                  5, // number of points
                                  3, // number of intervals
                                  2.0, // sigma_t
                                  0.8, // sigma_s
                                  1.1, // nu_sigma_f
                                  1.0, // internal source
                                  0.0, // boundary source
                                  1.0, // alpha
                                  2.0, // length
                                  1e-4); // tolerance

        // Test 2D steady state with boundary source
        cout << description << "steady state with boundary source, source iteration" << endl;
        checksum += test_infinite(true, // mls basis
                                  true, // mls weight
                                  "wendland11", // basis type
                                  "wendland11", // weight type
                                  weight_options,
                                  "source_iteration",
                                  2, // dimension
                                  2, // angular rule
                                  5, // number of points
                                  3, // number of intervals
                                  2.0, // sigma_t
                                  0.8, // sigma_s
                                  1.1, // nu_sigma_f
                                  1.0, // internal source
                                  1.0 / (2. * M_PI * (2.0 - 0.8 - 1.1)), // boundary source
                                  0.0, // alpha
                                  2.0, // length
                                  1e-4); // tolerance
    }

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
