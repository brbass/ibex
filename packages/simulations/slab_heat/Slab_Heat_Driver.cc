#include <cmath>
#include <iostream>
#include <fstream>
#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

#include "Cartesian_Plane.hh"
#include "Heat_Transfer_Integration.hh"
#include "Heat_Transfer_Factory.hh"
#include "Heat_Transfer_Solve.hh"
#include "Heat_Transfer_Solution.hh"
#include "Quadrature_Rule.hh"
#include "Slab_Heat_Data.hh"
#include "Solid_Geometry.hh"
#include "Timer.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;

void output_error(XML_Node output_node,
                  shared_ptr<Weak_Spatial_Discretization> spatial,
                  shared_ptr<Slab_Heat_Data> data,
                  shared_ptr<Heat_Transfer_Solution> solution)
{
    // Get number of integration ordinates
    int dimension = spatial->dimension();
    vector<int> dimensional_points(dimension);
    for (int i = 0; i < dimension; ++i)
    {
        int points_per_cell = 4;
        dimensional_points[i] = points_per_cell * spatial->options()->dimensional_cells[i];
    }
    
    // Get integration quadrature
    vector<vector<double> > limits = spatial->options()->limits;
    vector<vector<double> > ordinates;
    vector<double> weights;
    bool quad_available
        = Quadrature_Rule::cartesian_nd(dimension,
                                        Quadrature_Rule::Quadrature_Type::GAUSS_LEGENDRE,
                                        dimensional_points,
                                        limits,
                                        ordinates,
                                        weights);
    int number_of_ordinates = weights.size();

    // Get L_2 error
    double l2_num = 0;
    double l2_den = 0;
    for (int i = 0; i < number_of_ordinates; ++i)
    {
        vector<double> position = ordinates[i];
        double analytic = data->get_solution(position);
        double numeric = solution->solution(position);
        double diff = analytic - numeric;
        l2_num += diff * diff;
        l2_den += analytic * analytic;
    }
    double l2 = sqrt(l2_num / l2_den);
    output_node.set_child_value(l2, "l2_error");

    // Get collocation results
    int number_of_points = spatial->number_of_points();
    vector<double> values(number_of_points);
    vector<double> error(number_of_points);
    vector<double> benchmark(number_of_points);
    for (int i = 0; i < number_of_points; ++i)
    {
        vector<double> position = spatial->weight(i)->position();
        values[i] = solution->solution(position);
        benchmark[i] = data->get_solution(position);
    }
    output_node.set_child_vector(values, "solution");
    output_node.set_child_vector(benchmark, "benchmark");
}

void run_test(XML_Node input_node,
              XML_Node output_node)
{
    // Get dimension
    int dimension = input_node.get_child_value<int>("dimension");
    
    // Get solid geometry
    vector<vector<double> > limits = input_node.get_child_matrix<double>("limits",
                                                                         dimension,
                                                                         2); // surfaces per dimension
    shared_ptr<Solid_Geometry> solid;
    vector<shared_ptr<Cartesian_Plane> > surfaces;
    Heat_Transfer_Factory factory;
    factory.get_solid(dimension,
                      limits,
                      solid,
                      surfaces);

    // Get weak spatial discretization
    Weak_Spatial_Discretization_Parser spatial_parser(solid,
                                                      surfaces);
    shared_ptr<Weak_Spatial_Discretization> spatial
        = spatial_parser.get_weak_discretization(input_node.get_child("spatial_discretization"));

    // Get heat transfer integration options
    XML_Node heat_node = input_node.get_child("heat");
    string geometry = heat_node.get_attribute<string>("geometry");
    shared_ptr<Heat_Transfer_Integration_Options> integration_options
        = make_shared<Heat_Transfer_Integration_Options>();
    if (geometry == "cylindrical")
    {
        integration_options->geometry = Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_1D;
        AssertMsg(false, "analytic solution not implemented for cylindrical");
    }
    else if (geometry == "cartesian")
    {
        integration_options->geometry = Heat_Transfer_Integration_Options::Geometry::CARTESIAN;
    }
    else
    {
        AssertMsg(false, "geometry (" + geometry + ") not found");
    }
    
    // Get heat transfer data
    vector<double> conduction
        = heat_node.get_child_vector<double>("conduction",
                                             2);
    vector<double> convection
        = heat_node.get_child_vector<double>("convection",
                                             2);
    vector<double> source
        = heat_node.get_child_vector<double>("source",
                                             3);
    vector<double> temperature_inf
        = heat_node.get_child_vector<double>("temperature_inf",
                                             2);
    shared_ptr<Slab_Heat_Data> data
        = make_shared<Slab_Heat_Data>(integration_options,
                                      conduction,
                                      source,
                                      convection,
                                      temperature_inf,
                                      limits[0]);

    // Perform integration
    shared_ptr<Heat_Transfer_Integration> integration
        = make_shared<Heat_Transfer_Integration>(integration_options,
                                                 data,
                                                 spatial);

    // Initialize heat transfer solver
    shared_ptr<Heat_Transfer_Solve> solver
        = make_shared<Heat_Transfer_Solve>(integration,
                                           spatial);

    // Solve problem
    shared_ptr<Heat_Transfer_Solution> solution
        = solver->solve();

    // Output error and solution
    output_error(output_node.append_child("solution"),
                 spatial,
                 data,
                 solution);
    
    // Output other data
    solid->output(output_node.append_child("solid_geometry"));
    spatial->output(output_node.append_child("spatial_discretization"));
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

    // Get input and output nodes
    string input_filename = argv[1];
    XML_Document input_file(input_filename);
    XML_Node input_node = input_file.get_child("input");
    string output_filename = input_filename + ".out";
    XML_Document output_file;
    XML_Node output_node = output_file.append_child("output");
    run_test(input_node,
             output_node);
    output_file.save(output_filename);
    
    // Close MPI
    MPI_Finalize();

    return 0;
}
