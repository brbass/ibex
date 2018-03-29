#include <iostream>
#include <fstream>
#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source_Parser.hh"
#include "Cartesian_Plane.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Constructive_Solid_Geometry_Parser.hh"
#include "Cross_Section.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "LDFE_Quadrature.hh"
#include "Material.hh"
#include "Material_Parser.hh"
#include "Meshless_Function_Factory.hh"
#include "Quadrature_Rule.hh"
#include "Timer.hh"
#include "Transport_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;

shared_ptr<Weak_Spatial_Discretization>
get_spatial(XML_Node input_node)
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
    vector<shared_ptr<Material> > materials
        = material_parser.parse_from_xml(input_node.get_child("materials"));
    
    // Get boundary source
    Boundary_Source_Parser boundary_parser(angular,
                                           energy);
    vector<shared_ptr<Boundary_Source> > boundary_sources
        = boundary_parser.parse_from_xml(input_node.get_child("boundary_sources"));
    Assert(boundary_sources.size() == 1);
    
    // Get solid geometry
    Constructive_Solid_Geometry_Parser solid_parser(materials,
                                                    boundary_sources);
    shared_ptr<Constructive_Solid_Geometry> solid
        = solid_parser.parse_from_xml(input_node_.get_child("solid_geometry"));
    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces
        = solid->cartesian_boundary_surfaces();
    
    // Get boundary surfaces
    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces
        = solid->cartesian_boundary_surfaces();
    
    // Get spatial discretization
    Weak_Spatial_Discretization_Parser spatial_parser(solid,
                                                      boundary_surfaces);
    return spatial_parser.get_weak_discretization(input_node.get_child("spatial_discretization"));
}

// shared_ptr<VERA_Temperature> 
// run_heat(XML_Node input_node,
//          shared_ptr<VERA_Transport_Result> result,
//          shared_ptr<VERA_Temperature> weighting_temperature)
// {
//     // Get solid geometry
//     double length = 0.475;
//     shared_ptr<Solid_Geometry> solid;
//     vector<shared_ptr<Cartesian_Plane> > surfaces;
//     Heat_Transfer_Factory factory;
//     factory.get_solid(1, // dimension
//                       {{0, length}}, // limits
//                       solid,
//                       surfaces);

//     // Get weak spatial discretization
//     Weak_Spatial_Discretization_Parser spatial_parser(solid,
//                                                       surfaces);
//     shared_ptr<Weak_Spatial_Discretization> spatial
//         = spatial_parser.get_weak_discretization(input_node.get_child("spatial_discretization"));
    
//     // Get heat transfer data
//     shared_ptr<VERA_Heat_Data> data
//         = make_shared<VERA_Heat_Data>(result,
//                                       weighting_temperature);

//     // Get heat transfer integration
//     shared_ptr<Heat_Transfer_Integration_Options> integration_options
//         = make_shared<Heat_Transfer_Integration_Options>();
//     integration_options->geometry = Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_1D;
//     shared_ptr<Heat_Transfer_Integration> integration
//         = make_shared<Heat_Transfer_Integration>(integration_options,
//                                                  data,
//                                                  spatial);
    
//     // Get heat transfer solver
//     shared_ptr<Heat_Transfer_Solve> solver
//         = make_shared<Heat_Transfer_Solve>(integration,
//                                            spatial);
//     shared_ptr<Heat_Transfer_Solution> solution
//         = solver->solve();
    
//     return make_shared<VERA_Temperature>([solution](vector<double> const & position) -> double
//                                          {
//                                              double radius = sqrt(position[0] * position[0] + position[1] * position[1]);
//                                              if (radius > 0.475)
//                                              {
//                                                  return 600.0;
//                                              }
//                                              else
//                                              {
//                                                  return solution->solution({radius});
//                                              }
//                                          });
// }

// void output_temperature(shared_ptr<VERA_Temperature> temperature,
//                         XML_Node output_node)
// {
//     double dr = 0.001 - 1e-12;
//     vector<double> positions;
//     vector<double> values;
//     double maxr = 0.475 + 0.5 * dr;
//     int num_values = floor(maxr/dr) + 1;
//     for (int i = 0; i < num_values; ++i)
//     {
//         double r = i * dr;
//         vector<double> position = {r, 0};
//         positions.push_back(r);
//         values.push_back((*temperature)(position));
//     }
//     // for (double r = 0; r < maxr; r += dr)
//     // {
//     //     vector<double> position = {r, 0};
//     //     positions.push_back(r);
//     //     values.push_back((*temperature)(position));
//     // }

//     output_node.set_child_vector(positions, "points");
//     output_node.set_child_vector(values, "values");
// }

// void run_test(XML_Node input_node,
//               XML_Node output_node)
// {
//     // Get initial temperature
//     shared_ptr<VERA_Temperature> temperature
//         = make_shared<VERA_Temperature>([](vector<double> const &){return 600;});

//     // Get result pointer
//     shared_ptr<VERA_Transport_Result> result;

//     // Get total desired power
//     double pincell_power
//         = input_node.get_child("heat").get_child_value<double>("pincell_power");
    
//     Timer timer;
//     timer.start();
//     int num_iters = 4;
//     vector<double> eigenvalue_history;
//     for (int i = 0; i < num_iters; ++i)
//     {
//         // Run transport calculation
//         cout << "start transport calculation " << i << endl;
//         result
//             = run_transport(pincell_power,
//                             input_node.get_child("transport"),
//                             temperature);
//         eigenvalue_history.push_back(result->result()->k_eigenvalue);
//         cout << "end transport calculation " << i << endl;
        
//         // Run heat transfer calculation
//         cout << "start heat transfer calculation " << i << endl;
//         shared_ptr<VERA_Temperature> old_temperature = temperature;
//         temperature
//             = run_heat(input_node.get_child("heat"),
//                        result,
//                        old_temperature);
//         cout << "end heat transfer calculation " << i << endl;
//     }
    
//     // Output data
//     timer.stop();
//     output_temperature(temperature,
//                        output_node.append_child("temperature"));
//     output_node.append_child("timing").set_child_value(timer.time(), "total");
//     output_node.set_child_value(pincell_power, "pincell_power");
//     output_node.set_child_vector(eigenvalue_history, "eigenvalue_by_iteration");
//     result->output_data(output_node);
// }

void output_values(XML_Node input_node,
                   XML_Node output_node)
{
    // Get spatial discretization
    shared_ptr<Weak_Spatial_Discretization> spatial
        = get_spatial(input_node);
    int number_of_points = spatial->number_of_points();
    
    // Get grid of points for evaluation
    int dimension = spatial->dimension();
    XML_Node points_node = input_node.get_child("output_points");
    vector<vector<double> > limits = points_node.get_child_matrix<double>("limits",
                                                                          dimension,
                                                                          2); // surfaces
    vector<int> dimensional_points = points_node.get_child_vector<int>("dimensional_points",
                                                                       dimension);
    Meshless_Function_Factory factory;
    int number_of_evaluation_points;
    vector<vector<double> > evaluation_points;
    factory.get_cartesian_points(dimension,
                                 dimensional_points,
                                 limits,
                                 number_of_evaluation_points,
                                 evaluation_points);

    // Get values of basis functions at points
    vector<double> values(number_of_points * number_of_evaluation_points, 0);
    for (int i = 0; i < number_of_evaluation_points; ++i)
    {
        // Get nearest weight function
        vector<double> position = evaluation_points[i];
        int index = spatial->nearest_point(position);
        shared_ptr<Weight_Function> weight = spatial->weight(i);

        // Get values of all basis functions at this point
        int number_of_basis_functions = weight->number_of_basis_functions();
        vector<int> const basis_indices = weight->basis_function_indices();
        vector<double> basis_vals(number_of_basis_functions);
        for (int j = 0; j < number_of_basis_functions; ++j)
        {
            basis_vals[j] = weight->basis_function(j)->function()->base_function()->value(position);
        }

        // Normalize if applicable
        if (basis_depends_on_neighbors_)
        {
            vector<vector<double> > center_positions(number_of_basis_functions);
            for (int j = 0; j < number_of_basis_functions; ++j)
            {
                center_positions[j] = weight->basis_function(j)->position();
            }
            shared_ptr<Meshless_Normalization> norm
                = weight->basis_function(0)->function()->normalization();
            norm->get_values(position,
                             center_positions,
                             basis_vals,
                             basis_vals);
        }

        // Put the basis function values into the global matrix
        for (int j = 0; j < number_of_basis_functions; ++j)
        {
            int basis_index = basis_indices[j];
            int k = i + number_of_evaluation_points * basis_index;
            values[k] = basis_vals[j];
        }
    }

    // Get flattened points
    vector<double> flattened_evaluation_points(dimension * number_of_evaluation_points);
    for (int i = 0; i < number_of_evaluation_points; ++i)
    {
        for (int d = 0; d < dimension; ++d)
        {
            flattened_evaluation_points[d + dimension * i] = evaluation_points[i][d];
        }
    }
    vector<double> flattened_points(dimension * number_of_points);
    for (int i = 0; i < number_of_points; ++i)
    {
        for (int d = 0; d < dimension; ++d)
        {
            flattened_points[d + dimension * i] = points[i][d];
        }
    }

    output_node.set_child_vector(flattened_evaluation_points, "evaluation_points");
    output_node.set_child_vector(values, "evaluation_values");
    output_node.set_child_vector(flattened_points, "points");
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
    output_values(input_node,
                  output_node);
    output_file.save(output_filename);
    
    // Close MPI
    MPI_Finalize();

    return 0;
}
