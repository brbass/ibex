#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source_Parser.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Constructive_Solid_Geometry_Parser.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Material.hh"
#include "Material_Parser.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;

// Get number of angles that reach the first region in the solid geometry for given points
void output_values(XML_Node input_node,
                   XML_Node output_node)
{
    // Get data
    Energy_Discretization_Parser energy_parser;
    shared_ptr<Energy_Discretization> energy = 
        energy_parser.parse_from_xml(input_node.get_child("energy_discretization"));
    
    Angular_Discretization_Parser angular_parser;
    shared_ptr<Angular_Discretization> angular = 
        angular_parser.parse_from_xml(input_node.get_child("angular_discretization"));
    
    Material_Parser material_parser(angular,
                                    energy);
    vector<shared_ptr<Material> > materials
        = material_parser.parse_from_xml(input_node.get_child("materials"));
    
    Boundary_Source_Parser boundary_parser(angular,
                                           energy);
    vector<shared_ptr<Boundary_Source> > boundary_sources
        = boundary_parser.parse_from_xml(input_node.get_child("boundary_sources"));
    
    Constructive_Solid_Geometry_Parser solid_parser(materials,
                                                    boundary_sources);
    shared_ptr<Constructive_Solid_Geometry> solid
        = solid_parser.parse_from_xml(input_node.get_child("solid_geometry"));

    int dimension = angular->dimension();
    int number_of_ordinates = angular->number_of_ordinates();
    XML_Node eval_node = input_node.get_child("evaluation");
    int number_of_evaluation_points = eval_node.get_child_value<int>("number");
    vector<vector<double> > eval_points = eval_node.get_child_matrix<double>("points",
                                                                             number_of_evaluation_points,
                                                                             dimension);
    
    // Get number of vectors from the source point that hit the first region
    vector<int> num_hits(number_of_evaluation_points, 0);
    for (int i = 0; i < number_of_evaluation_points; ++i)
    {
        vector<double> position = eval_points[i];
        
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            vector<double> direction = angular->direction(o);

            int final_region;
            double distance;
            vector<double> final_position;
            int surface = solid->next_intersection(position,
                                                   direction,
                                                   final_region,
                                                   distance,
                                                   final_position);

            if (final_region == 0)
            {
                num_hits[i] += 1;
            }
        }
    }
    
    output_node.set_child_value(dimension, "dimension");
    output_node.set_child_value(number_of_ordinates, "number_of_ordinates");
    output_node.set_child_value(number_of_evaluation_points, "number_of_evaluation_points");
    output_node.set_child_matrix(eval_points, "evaluation_points");
    output_node.set_child_vector(num_hits, "number_of_hits");
}

int main(int argc, char **argv)
{
    string input_filename = argv[1];
    XML_Document input_file(input_filename);
    XML_Node input_node = input_file.get_child("input");
    string output_filename = input_filename + ".out";
    XML_Document output_file;
    XML_Node output_node = output_file.append_child("output");
    output_values(input_node,
                  output_node);
    output_file.save(output_filename);
    
    return 0;
}
