#include <cmath>
#include <iostream>
#include <fstream>
#include <memory>
#include <mpi.h>
#if defined(ENABLE_OPENMP)
    #include <omp.h>
#else
    inline void omp_set_num_threads(int i) {return;}
    inline void omp_get_num_threads() {return 1;}
#endif
#include <string>
#include <vector>

#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source_Parser.hh"
#include "Cartesian_Plane.hh"
#include "Cross_Section.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Heat_Transfer_Integration.hh"
#include "Heat_Transfer_Factory.hh"
#include "Heat_Transfer_Solve.hh"
#include "Heat_Transfer_Solution.hh"
#include "Krylov_Eigenvalue.hh"
#include "LDFE_Quadrature.hh"
#include "Material.hh"
#include "Material_Parser.hh"
#include "Meshless_Sweep.hh"
#include "Meshless_Sweep_Parser.hh"
#include "Quadrature_Rule.hh"
#include "Random_Number_Generator.hh"
#include "Solver.hh"
#include "Solver_Parser.hh"
#include "Timer.hh"
#include "Transport_Discretization.hh"
#include "VERA_Heat_Data.hh"
#include "VERA_Solid_Geometry.hh"
#include "VERA_Transport_Result.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;

shared_ptr<VERA_Transport_Result>
run_transport(bool include_crack,
              int heat_dimension,
              double pincell_power,
              XML_Node input_node,
              shared_ptr<VERA_Temperature> temperature,
              int num_threads)
{
    // Get energy discretization
    Energy_Discretization_Parser energy_parser;
    shared_ptr<Energy_Discretization> energy
        = energy_parser.parse_from_xml(input_node.get_child("energy_discretization"));

    // Get angular discretization
    Angular_Discretization_Parser angular_parser;
    shared_ptr<Angular_Discretization> angular
        = angular_parser.parse_from_xml(input_node.get_child("angular_discretization"));

    // Get materials
    Material_Parser material_parser(angular,
                                    energy);
    bool include_ifba = input_node.get_child("materials").get_attribute<bool>("include_ifba");
    vector<shared_ptr<Material> > materials
        = material_parser.parse_from_xml(input_node.get_child("materials"));
    
    // Get boundary source
    Boundary_Source_Parser boundary_parser(angular,
                                           energy);
    vector<shared_ptr<Boundary_Source> > boundary_sources
        = boundary_parser.parse_from_xml(input_node.get_child("boundary_sources"));
    Assert(boundary_sources.size() == 1);
    
    // Get solid geometry
    shared_ptr<VERA_Solid_Geometry> solid
        = make_shared<VERA_Solid_Geometry>(include_ifba,
                                           include_crack,
                                           temperature,
                                           angular,
                                           energy,
                                           materials,
                                           boundary_sources[0]);
    
    // Get boundary surfaces
    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces
        = solid->cartesian_boundary_surfaces();
    
    // Get spatial discretization
    Weak_Spatial_Discretization_Parser spatial_parser(solid,
                                                      boundary_surfaces);
    // omp_set_num_threads(num_threads);
    shared_ptr<Weak_Spatial_Discretization> spatial
        = spatial_parser.get_weak_discretization(input_node.get_child("spatial_discretization"));
    // omp_set_num_threads(1);
    
    // Get transport discretization
    shared_ptr<Transport_Discretization> transport
        = make_shared<Transport_Discretization>(spatial,
                                                angular,
                                                energy);

    // Get sweep
    Meshless_Sweep_Parser sweep_parser(spatial,
                                   angular,
                                   energy,
                                   transport);
    // omp_set_num_threads(num_threads);
    shared_ptr<Meshless_Sweep> sweep
        = sweep_parser.get_meshless_sweep(input_node.get_child("transport"));
    // omp_set_num_threads(1);
    
    // Get solver
    Solver_Parser solver_parser(spatial,
                                angular,
                                energy,
                                transport);
    shared_ptr<Solver> solver
        = solver_parser.get_krylov_eigenvalue(input_node.get_child("solver"),
                                              sweep);

    // omp_set_num_threads(num_threads);
    solver->solve();
    // omp_set_num_threads(1);

    double fuel_radius;
    if (include_crack)
    {
        fuel_radius = 0.4135;
    }
    else
    {
        fuel_radius = 0.4096;
    }
    
    return make_shared<VERA_Transport_Result>(heat_dimension,
                                              fuel_radius,
                                              pincell_power,
                                              solid,
                                              angular,
                                              energy,
                                              spatial,
                                              solver,
                                              solver->result());
}

shared_ptr<VERA_Temperature> 
run_heat(bool include_crack,
         int heat_dimension,
         XML_Node input_node,
         shared_ptr<VERA_Transport_Result> result,
         shared_ptr<VERA_Temperature> weighting_temperature)
{
    if (include_crack)
    {
        Assert(heat_dimension == 2);
    }
    
    // Get solid geometry
    double length = 0.475;
    shared_ptr<Solid_Geometry> solid;
    vector<shared_ptr<Cartesian_Plane> > surfaces;
    Heat_Transfer_Factory factory;
    switch (heat_dimension)
    {
    case 1:
        factory.get_solid(heat_dimension,
                          {{0, length}}, // limits
                          solid,
                          surfaces);
        break;
    case 2:
        factory.get_solid(heat_dimension,
                          {{-length, length}, {-length, length}}, // limits
                          solid,
                          surfaces);
        break;
    default:
        AssertMsg(false, "dimension not found");
    }

    // Get weak spatial discretization
    Weak_Spatial_Discretization_Parser spatial_parser(solid,
                                                      surfaces);
    shared_ptr<Weak_Spatial_Discretization> spatial
        = spatial_parser.get_weak_discretization(input_node.get_child("spatial_discretization"));
    
    // Get heat transfer data
    shared_ptr<VERA_Heat_Data> data
        = make_shared<VERA_Heat_Data>(include_crack,
                                      heat_dimension,
                                      result,
                                      weighting_temperature);

    // Get heat transfer integration
    shared_ptr<Heat_Transfer_Integration_Options> integration_options
        = make_shared<Heat_Transfer_Integration_Options>();
    switch (heat_dimension)
    {
    case 1:
        integration_options->geometry = Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_1D;
        break;
    case 2:
        integration_options->geometry = Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_2D;
        break;
    default:
        AssertMsg(false, "dimension not found");
    }
    
    shared_ptr<Heat_Transfer_Integration> integration
        = make_shared<Heat_Transfer_Integration>(integration_options,
                                                 data,
                                                 spatial);
    
    // Get heat transfer solver
    shared_ptr<Heat_Transfer_Solve> solver
        = make_shared<Heat_Transfer_Solve>(integration,
                                           spatial);
    shared_ptr<Heat_Transfer_Solution> solution
        = solver->solve();

    switch (heat_dimension)
    {
    case 1:
        
        return make_shared<VERA_Temperature>([solution](vector<double> const &position) -> double
            {
                double radius2 = position[0] * position[0] + position[1] * position[1];
                if (radius2 > 0.475 * 0.475 + 1e-12)
                {
                    return 600.0;
                }
                else
                {
                    return solution->solution({sqrt(radius2)});
                }
            });
    case 2:
        return make_shared<VERA_Temperature>([solution](vector<double> const &position) -> double
            {
                double radius2 = position[0] * position[0] + position[1] * position[1];
                if (radius2 > 0.475 * 0.475 + 1e-12)
                {
                    return 600.0;
                }
                else
                {
                    return solution->solution(position);
                }
            });
    default:
        AssertMsg(false, "dimension not found");
    }
}

void output_temperature(int heat_dimension,
                        shared_ptr<VERA_Temperature> temperature,
                        XML_Node output_node)
{
    vector<double> positions;
    vector<double> values;
    
    switch (heat_dimension)
    {
    case 1:
    {
        double dr = 0.001 - 1e-12;
        double maxr = 0.475 + 0.5 * dr;
        int num_values = floor(maxr/dr) + 1;
        for (int i = 0; i < num_values; ++i)
        {
            double r = i * dr;
            vector<double> position = {r, 0};
            positions.push_back(r);
            values.push_back((*temperature)(position));
        }
        break;
    }
    case 2:
    {
        double dr = 0.01 - 1e-11;
        double maxr = 0.475 + 0.5 * dr;
        int num_values = floor(maxr/dr) + 1;
        Random_Number_Generator<double> rng(0,
                                            2 * M_PI,
                                            492); // seed
        for (int i = 0; i < num_values; ++i)
        {
            double r = i * dr;

            double circ = 2 * M_PI * r;
            int numt = int(ceil(circ / dr));
            double start_t = rng.scalar();
            double dt = 2 * M_PI / static_cast<double>(numt);
            
            for (int i = 0; i < numt; ++i)
            {
                double t = start_t + i * dt;
                vector<double> position = {r * cos(t), r * sin(t)};
                for (int d = 0; d < heat_dimension; ++d)
                {
                    positions.push_back(position[d]);
                }
                values.push_back((*temperature)(position));
            }
        }
        break;
    }
    default:
        AssertMsg(false, "dimension not found");
    }
    
    output_node.set_child_vector(positions, "points");
    output_node.set_child_vector(values, "values");
}

void run_test(XML_Node input_node,
              XML_Node output_node,
              int num_threads)
{
    // Get initial temperature
    shared_ptr<VERA_Temperature> temperature
        = make_shared<VERA_Temperature>([](vector<double> const &){return 600;});

    // Get result pointer
    shared_ptr<VERA_Transport_Result> result;

    // Get total desired power
    bool include_crack
        = input_node.get_child("heat").get_attribute<bool>("include_crack");
    int heat_dimension
        = input_node.get_child("heat").get_child_value<int>("heat_dimension");
    double pincell_power
        = input_node.get_child("heat").get_child_value<double>("pincell_power");
    int number_of_iterations
        = input_node.get_child("heat").get_child_value<int>("number_of_iterations");
    Timer timer;
    timer.start();
    vector<double> eigenvalue_history;
    for (int i = 0; i < number_of_iterations; ++i)
    {
        // Run transport calculation
        cout << "start transport calculation " << i << endl;
        result
            = run_transport(include_crack,
                            heat_dimension,
                            pincell_power,
                            input_node.get_child("transport"),
                            temperature,
                            num_threads);
        eigenvalue_history.push_back(result->result()->k_eigenvalue);
        cout << "end transport calculation " << i << endl;
        
        // Run heat transfer calculation
        cout << "start heat transfer calculation " << i << endl;
        shared_ptr<VERA_Temperature> old_temperature = temperature;
        temperature
            = run_heat(include_crack,
                       heat_dimension,
                       input_node.get_child("heat"),
                       result,
                       old_temperature);
        cout << "end heat transfer calculation " << i << endl;
    }
    
    // Output data
    timer.stop();
    output_temperature(heat_dimension,
                       temperature,
                       output_node.append_child("temperature"));
    output_node.append_child("timing").set_child_value(timer.time(), "total");
    output_node.set_child_value(pincell_power, "pincell_power");
    output_node.set_child_vector(eigenvalue_history, "eigenvalue_by_iteration");
    result->output_data(output_node);
}

int main(int argc, char **argv)
{
    // Initialize MPI
    MPI_Init(&argc, &argv);

    if (argc != 2)
    {
        cerr << "need input file" << endl;
        return 1;
    }
    
    // Get base XML node
    string input_filename = argv[1];
    string output_filename = input_filename + ".out";
    XML_Document input_file(input_filename);
    XML_Node input_node = input_file.get_child("input");
    XML_Document output_file;
    XML_Node output_node = output_file.append_child("output");

    // Set number of threads for transport
    int number_of_threads = input_node.get_attribute<int>("number_of_threads",
                                                          1);
    omp_set_num_threads(number_of_threads);
    
    // Run problem
    run_test(input_node,
             output_node,
             number_of_threads);

    // Output data
    output_file.save(output_filename);
    
    // Close MPI
    MPI_Finalize();

    return 0;
}
