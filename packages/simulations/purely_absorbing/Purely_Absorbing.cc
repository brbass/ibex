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
#include "Meshless_Sweep.hh"
#include "Meshless_Sweep_Parser.hh"
#include "SUPG_Internal_Source_Operator.hh"
#include "SUPG_Moment_To_Discrete.hh"
#include "Transport_Discretization.hh"
#include "Vector_Operator_Functions.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;

// Get number of angles that reach the first region in the solid geometry for given points
void run_problem(XML_Node input_node,
                 XML_Node output_node)
{
    // Energy discretization
    Energy_Discretization_Parser energy_parser;
    shared_ptr<Energy_Discretization> energy = 
        energy_parser.parse_from_xml(input_node.get_child("energy_discretization"));
    Assert(energy->number_of_groups() == 1);

    // Angular discretization
    Angular_Discretization_Parser angular_parser;
    shared_ptr<Angular_Discretization> angular = 
        angular_parser.parse_from_xml(input_node.get_child("angular_discretization"));

    // Material
    Material_Parser material_parser(angular,
                                    energy);
    vector<shared_ptr<Material> > materials
        = material_parser.parse_from_xml(input_node.get_child("materials"));
    Assert(materials.size() == 1);

    // Boundary source
    Boundary_Source_Parser boundary_parser(angular,
                                           energy);
    vector<shared_ptr<Boundary_Source> > boundary_sources
        = boundary_parser.parse_from_xml(input_node.get_child("boundary_sources"));

    // Constructive solid geometry
    Constructive_Solid_Geometry_Parser solid_parser(materials,
                                                    boundary_sources);
    shared_ptr<Constructive_Solid_Geometry> solid
        = solid_parser.parse_from_xml(input_node.get_child("solid_geometry"));

    // Get boundary surfaces
    Assert(solid->cartesian_boundaries());
    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces
        = solid->cartesian_boundary_surfaces();

    // Get spatial discretization
    Weak_Spatial_Discretization_Parser spatial_parser(solid,
                                                      boundary_surfaces);
    shared_ptr<Weak_Spatial_Discretization>spatial
        = spatial_parser.get_weak_discretization(input_node.get_child("spatial_discretization"));
    Assert(spatial->options()->include_supg);
    Assert(!(spatial->has_reflection()));
    
    // Get transport discretization
    shared_ptr<Transport_Discretization> transport
        = make_shared<Transport_Discretization>(spatial,
                                                angular,
                                                energy);
    
    // Get sweep operator
    Meshless_Sweep_Parser sweep_parser(spatial,
                                       angular,
                                       energy,
                                       transport);
    shared_ptr<Meshless_Sweep> Linv
        = sweep_parser.get_meshless_sweep(input_node.get_child("transport"));
    Linv->set_include_boundary_source(true);
    
    // Get internal source operator
    shared_ptr<SUPG_Internal_Source_Operator> Q
        = make_shared<SUPG_Internal_Source_Operator>(spatial,
                                                     angular,
                                                     energy);
    
    // Get moment to discrete operator
    shared_ptr<SUPG_Moment_To_Discrete> M
        = make_shared<SUPG_Moment_To_Discrete>(spatial,
                                               angular,
                                               energy,
                                               false); // include double dimensional moments
    
    // Get combined operator
    shared_ptr<Vector_Operator> combined = Linv * M * Q;
    
    // Get solution vector
    vector<double> solution(Q->column_size());

    // Get solution
    (*combined)(solution);

    cout << "hi" << endl;
}

int main(int argc, char **argv)
{
    string input_filename = argv[1];
    XML_Document input_file(input_filename);
    XML_Node input_node = input_file.get_child("input");
    string output_filename = input_filename + ".out";
    XML_Document output_file;
    XML_Node output_node = output_file.append_child("output");
    run_problem(input_node,
                output_node);
    output_file.save(output_filename);
    
    return 0;
}
