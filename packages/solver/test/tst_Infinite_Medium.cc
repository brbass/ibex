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
#include "Constructive_Solid_Geometry.hh"
#include "Constructive_Solid_Geometry_Parser.hh"
#include "Cross_Section.hh"
#include "Discrete_Value_Operator.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
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

using namespace std;

void get_one_region(bool basis_mls,
                    bool weight_mls,
                    string basis_type,
                    string weight_type,
                    Weight_Function::Options weight_options,
                    bool krylov,
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
                                                    weight_options);
    
    // Get transport discretization
    transport
        = make_shared<Transport_Discretization>(spatial,
                                                angular,
                                                energy);
    
    // Get weak RBF sweep
    Weak_RBF_Sweep::Options options;
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
    
    if (krylov)
    {
        solver
            = solver_factory.get_krylov_steady_state(sweeper,
                                                     convergence);
    }
    else
    {
        solver
            = solver_factory.get_source_iteration(sweeper,
                                                  convergence);
    }
}

int test_infinite(bool basis_mls,
                  bool weight_mls,
                  string basis_type,
                  string weight_type,
                  Weight_Function::Options weight_options,
                  bool krylov,
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
                  double length)
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
                   krylov,
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
    
    // Print
    int phi_size = transport->phi_size();
    int num_values = result->phi.size();
    int w = 16;
    cout << setw(w) << "value" << setw(w) << "weighted" << endl;
    for (int i = 0; i < phi_size; ++i)
    {
        for (int j = 0; j < num_values; ++j)
        {
            cout << setw(w) << result->phi[j][i];
        }
        cout << endl;
    }

    return checksum;
}


int main(int argc, char **argv)
{
    int checksum = 0;
    
    MPI_Init(&argc, &argv);
    
    {
        bool krylov = true;
        int dimension = 1;
        int num_dimensional_points = 5;
        int angular_rule = dimension == 1 ? 128 : 3;
        double norm = dimension == 1 ? 2 : 2 * M_PI;
        double radius_num_intervals = 3.0;
        double sigma_t = 2.0;
        double sigma_s = 1.0;
        double chi_nu_sigma_f = 0.5;
        double internal_source = 1.0;
        double boundary_source = 1 / (norm * (sigma_t - sigma_s - chi_nu_sigma_f));
        double alpha = 0.0;
        double length = 10;
        bool basis_mls = true;
        bool weight_mls = true;
        string basis_type = "wendland11";
        string weight_type = "wendland11";
        Weight_Function::Options weight_options;
        weight_options.integration_ordinates = 32;
        weight_options.tau_const = 1.0;
        weight_options.tau_scaling = Weight_Function::Options::Tau_Scaling::NONE;
        // weight_options.output = Weight_Function::Options::Output::STANDARD;
        // checksum += test_infinite(basis_mls,
        //                           weight_mls,
        //                           basis_type,
        //                           weight_type,
        //                           weight_options,
        //                           krylov,
        //                           dimension,
        //                           angular_rule,
        //                           num_dimensional_points,
        //                           radius_num_intervals,
        //                           sigma_t,
        //                           sigma_s,
        //                           chi_nu_sigma_f,
        //                           internal_source,
        //                           boundary_source,
        //                           alpha,
        //                           length);
        weight_options.output = Weight_Function::Options::Output::SUPG;
        checksum += test_infinite(basis_mls,
                                  weight_mls,
                                  basis_type,
                                  weight_type,
                                  weight_options,
                                  krylov,
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
                                  length);
    }
    
    MPI_Finalize();
    
    return checksum;
}
