#include <iostream>
#include <fstream>
#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

#include "Angular_Discretization_Parser.hh"
#include "Basis_Function.hh"
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
#include "Meshless_Function.hh"
#include "Meshless_Function_Factory.hh"
#include "Meshless_Normalization.hh"
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
    
    // Get solid geometry
    Constructive_Solid_Geometry_Parser solid_parser(materials,
                                                    boundary_sources);
    shared_ptr<Constructive_Solid_Geometry> solid
        = solid_parser.parse_from_xml(input_node.get_child("solid_geometry"));
    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces
        = solid->cartesian_boundary_surfaces();
    
    // Get spatial discretization
    Weak_Spatial_Discretization_Parser spatial_parser(solid,
                                                      boundary_surfaces);
    return spatial_parser.get_weak_discretization(input_node.get_child("spatial_discretization"));
}

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
        shared_ptr<Weight_Function> weight = spatial->weight(index);

        // Get values of all basis functions at this point
        int number_of_basis_functions = weight->number_of_basis_functions();
        vector<int> const basis_indices = weight->basis_function_indices();
        vector<double> basis_vals(number_of_basis_functions);
        for (int j = 0; j < number_of_basis_functions; ++j)
        {
            basis_vals[j] = weight->basis_function(j)->function()->base_function()->value(position);
        }

        // Normalize if applicable
        if (weight->basis_function(0)->function()->depends_on_neighbors())
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
            flattened_points[d + dimension * i] = spatial->basis(i)->position()[d];
        }
    }

    // Get radii
    vector<double> radii(number_of_points);
    for (int i = 0; i < number_of_points; ++i)
    {
        radii[i] = spatial->basis(i)->radius();
    }

    // Output data
    output_node.set_child_value(number_of_points, "number_of_points");
    output_node.set_child_value(number_of_evaluation_points, "number_of_evaluation_points");
    output_node.set_child_value(dimension, "dimension");
    output_node.set_child_vector(points_node.get_child_vector<double>("limits",
                                                                      dimension * 2), "limits");
    output_node.set_child_vector(dimensional_points, "dimensional_points");
    output_node.set_child_vector(flattened_evaluation_points, "evaluation_points");
    output_node.set_child_vector(values, "evaluation_values");
    output_node.set_child_vector(flattened_points, "points");
    output_node.set_child_vector(radii, "radii");
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
