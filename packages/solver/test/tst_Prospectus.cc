#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <sstream>

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Factory.hh"
#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source.hh"
#include "Boundary_Source_Parser.hh"
#include "Cartesian_Plane.hh"
#include "Check.hh"
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
#include "Meshless_Function_Factory.hh"
#include "Region.hh"
#include "Solver_Factory.hh"
#include "Source_Iteration.hh"
#include "Timer.hh"
#include "Transport_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Factory.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "Weak_RBF_Sweep.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

namespace ce = Check_Equality;
using namespace std;

void get_pincell(bool basis_mls,
                 bool weight_mls,
                 string basis_type,
                 string weight_type,
                 string method,
                 Weight_Function::Options weight_options,
                 Weak_RBF_Sweep::Options sweep_options,
                 int dimension,
                 int angular_rule,
                 int num_dimensional_points,
                 double radius_num_intervals,
                 shared_ptr<Weak_Spatial_Discretization> &spatial,
                 shared_ptr<Angular_Discretization> &angular,
                 shared_ptr<Energy_Discretization> &energy,
                 shared_ptr<Solver> &solver)
{
    // Set constants
    double length = dimension == 1 ? 2.0 : 4.0;
    
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
    vector<shared_ptr<Material> > materials;
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
                                                 {1.0}); // internal source
    materials[1]
        = material_factory.get_standard_material(1, // index
                                                 {2.0}, // sigma_t
                                                 {1.9}, // sigma_s
                                                 {0.0}, // nu
                                                 {0.0}, // sigma_f
                                                 {0.0}, // chi
                                                 {0.0}); // internal source
    
    // Get boundary source
    vector<shared_ptr<Boundary_Source> > boundary_sources;
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
    
    // Get transport discretization
    shared_ptr<Transport_Discretization> transport
        = make_shared<Transport_Discretization>(spatial,
                                                angular,
                                                energy);
    
    // Get weak RBF sweep
    shared_ptr<Weak_RBF_Sweep> sweeper
        = make_shared<Weak_RBF_Sweep>(sweep_options,
                                      spatial,
                                      angular,
                                      energy,
                                      transport);

    // Get convergence method
    shared_ptr<Convergence_Measure> convergence
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

void get_one_region(bool basis_mls,
                    bool weight_mls,
                    string basis_type,
                    string weight_type,
                    string method,
                    Weight_Function::Options weight_options,
                    Weak_RBF_Sweep::Options sweep_options,
                    int dimension,
                    int angular_rule,
                    int num_dimensional_points,
                    double radius_num_intervals,
                    shared_ptr<Weak_Spatial_Discretization> &spatial,
                    shared_ptr<Angular_Discretization> &angular,
                    shared_ptr<Energy_Discretization> &energy,
                    shared_ptr<Solver> &solver)
{
    double length = dimension == 1 ? 2.0 : 4.0;
    
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
    vector<shared_ptr<Material> > materials;
    materials.resize(1);
    Material_Factory material_factory(angular,
                                      energy);
    double sigma_t = 1.0;
    double sigma_s = 0.84;
    double nu = 2.4;
    double sigma_f = 0.1;
    double chi = 1.;
    double internal_source = 1.0;
    materials[0]
        = material_factory.get_standard_material(0, // index
                                                 {sigma_t},
                                                 {sigma_s}, // sigma_s
                                                 {nu}, // nu
                                                 {sigma_f}, // sigma_f
                                                 {chi}, // chi
                                                 {internal_source});
    
    // Get boundary source
    vector<shared_ptr<Boundary_Source> > boundary_sources;
    Boundary_Source::Dependencies boundary_dependencies;
    if (dimension == 1)
    {
        boundary_sources.resize(2);
        boundary_sources[0]
            = make_shared<Boundary_Source>(0, // index
                                           boundary_dependencies,
                                           angular,
                                           energy,
                                           vector<double>(number_of_groups, 0.), // boundary source
                                           vector<double>(number_of_groups, 1.0)); // alpha
        boundary_sources[1]
            = make_shared<Boundary_Source>(0, // index
                                           boundary_dependencies,
                                           angular,
                                           energy,
                                           vector<double>(number_of_groups, 0.), // boundary source
                                           vector<double>(number_of_groups, 0.)); // alpha
    }
    else
    {
        boundary_sources.resize(1);
        boundary_sources[0]
            = make_shared<Boundary_Source>(0, // index
                                           boundary_dependencies,
                                           angular,
                                           energy,
                                           vector<double>(number_of_groups, 0.), // boundary source
                                           vector<double>(number_of_groups, 0.0)); // alpha
    }
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
    if (dimension == 1)
    {
        surfaces[0]->set_boundary_source(boundary_sources[0]);
        surfaces[1]->set_boundary_source(boundary_sources[1]);
    }
    else
    {
        for (shared_ptr<Surface> surface : surfaces)
        {
            surface->set_boundary_source(boundary_sources[0]);
        }
    }
    vector<Surface::Relation> surface_relations(2 * dimension,
                                                Surface::Relation::NEGATIVE);
    regions[0]
        = make_shared<Region>(0, // index
                              materials[0],
                              surface_relations,
                              surfaces);
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
    
    // Get transport discretization
    shared_ptr<Transport_Discretization> transport
        = make_shared<Transport_Discretization>(spatial,
                                                angular,
                                                energy);
    
    // Get weak RBF sweep
    Weak_RBF_Sweep::Options options;
    shared_ptr<Weak_RBF_Sweep> sweeper
        = make_shared<Weak_RBF_Sweep>(options,
                                      spatial,
                                      angular,
                                      energy,
                                      transport);

    // Get convergence method
    shared_ptr<Convergence_Measure> convergence
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

void output_results(shared_ptr<Solver::Result> result,
                    shared_ptr<Weak_Spatial_Discretization> spatial,
                    double setup_time,
                    double solve_time,
                    string output_path)
{
    // Get XML document
    XML_Document output_document;
    XML_Node output_node = output_document.append_child("output");

    // Set scalar quantities
    output_node.set_child_value(result->k_eigenvalue, "k_eigenvalue");
    output_node.set_child_value(result->total_iterations, "total_iterations");
    output_node.set_child_value(result->source_iterations, "source_iterations");
    output_node.set_child_value(spatial->dimension(), "dimension");
    output_node.set_child_value(setup_time, "setup_time");
    output_node.set_child_value(solve_time, "solve_time");
    output_node.set_child_value(setup_time + solve_time, "total_time");
    
    // Get grid of points to calculate scalar flux
    int dimension = spatial->dimension();
    vector<int> dimensional_points(dimension, 100);
    vector<vector<double> > limits = spatial->options().limits;
    Meshless_Function_Factory factory;
    int number_of_result_points;
    vector<vector<double> > points;
    factory.get_cartesian_points(dimension,
                                 dimensional_points,
                                 limits,
                                 number_of_result_points,
                                 points);

    // Get values for flux at each point
    vector<double> scalar_flux(number_of_result_points);
    for (int i = 0; i < number_of_result_points; ++i)
    {
        vector<double> position = points[i];
        scalar_flux[i] = spatial->expansion_value(position,
                                                  result->coefficients);
    }

    // Flatten point array
    vector<double> flat_points(dimension * number_of_result_points);
    for (int i = 0; i < number_of_result_points; ++i)
    {
        for (int d = 0; d < dimension; ++d)
        {
            flat_points[d + dimension * i] = points[i][d];
        }
    }
    
    // Store values
    output_node.set_child_vector(flat_points, "points");
    output_node.set_child_vector(scalar_flux, "flux");
    
    // Save XML document
    output_document.save(output_path);
}

void run_problem(int test_num,
                 bool basis_mls,
                 bool weight_mls,
                 string basis_type,
                 string weight_type,
                 string method,
                 Weight_Function::Options weight_options,
                 Weak_RBF_Sweep::Options sweep_options,
                 string output_path,
                 int dimension,
                 int angular_rule,
                 int num_dimensional_points,
                 double radius_num_intervals)
{
    // Initialize timing
    Timer timer;
    
    // Get problem
    timer.start();
    shared_ptr<Weak_Spatial_Discretization> spatial;
    shared_ptr<Angular_Discretization> angular;
    shared_ptr<Energy_Discretization> energy;
    shared_ptr<Solver> solver;
    if (test_num == 1)
    {
        get_one_region(basis_mls,
                       weight_mls,
                       basis_type,
                       weight_type,
                       method,
                       weight_options,
                       sweep_options,
                       dimension,
                       angular_rule,
                       num_dimensional_points,
                       radius_num_intervals,
                       spatial,
                       angular,
                       energy,
                       solver);
    }
    else if (test_num == 2)
    {
        get_pincell(basis_mls,
                    weight_mls,
                    basis_type,
                    weight_type,
                    method,
                    weight_options,
                    sweep_options,
                    dimension,
                    angular_rule,
                    num_dimensional_points,
                    radius_num_intervals,
                    spatial,
                    angular,
                    energy,
                    solver);
    }
    else
    {
        AssertMsg(false, "test number not found");
    }
    timer.stop();
    double setup_time = timer.time();

    // Solve problem and get result
    timer.start();
    solver->solve();
    timer.stop();
    double solve_time = timer.time();
    shared_ptr<Solver::Result> result
        = solver->result();
    
    // Output results
    output_results(result,
                   spatial,
                   setup_time,
                   solve_time,
                   output_path);
}

void run_test(string input_path)
{
    // Get XML document
    string output_path = input_path + ".out";
    XML_Document input_file(input_path);
    XML_Node input_node = input_file.get_child("input");
    
    // Get arguments
    int test_num = input_node.get_child_value<int>("test_num");
    bool basis_mls = input_node.get_child_value<bool>("basis_mls");
    bool weight_mls = input_node.get_child_value<bool>("weight_mls");
    string basis_type = input_node.get_child_value<string>("basis_type");
    string weight_type = input_node.get_child_value<string>("weight_type");
    string method = input_node.get_child_value<string>("solver");
    bool supg = input_node.get_child_value<bool>("supg");
    int dimension = input_node.get_child_value<int>("dimension");
    int angular_rule = input_node.get_child_value<int>("angular_rule");
    int num_dimensional_points = input_node.get_child_value<int>("num_dimensional_points");
    double radius_num_intervals = input_node.get_child_value<double>("radius_num_intervals");

    // Set options
    Weight_Function::Options weight_options;
    if (supg)
    {
        weight_options.output = Weight_Function::Options::Output::SUPG;
    }
    else
    {
        weight_options.output = Weight_Function::Options::Output::STANDARD;
    }
    if (dimension == 1)
    {
        weight_options.integration_ordinates = 32;
    }
    else
    {
        weight_options.integration_ordinates = 16;
    }
    
    Weak_RBF_Sweep::Options sweep_options;
    sweep_options.solver = Weak_RBF_Sweep::Options::Solver::AMESOS;

    // Run problem
    run_problem(test_num,
                basis_mls,
                weight_mls,
                basis_type,
                weight_type,
                method,
                weight_options,
                sweep_options,
                output_path,
                dimension,
                angular_rule,
                num_dimensional_points,
                radius_num_intervals);
}    

int main(int argc, char **argv)
{
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get input document path
    if (argc != 2)
    {
        cout << "args: input_path" << endl;
        return 1;
    }
    string input_path = argv[1];
    
    // Run test
    run_test(input_path);
    
    // Close MPI
    MPI_Finalize();
    
    return 0;
}
