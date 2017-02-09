#include <memory>
#include <vector>

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Parser.hh"
#include "Basis_Function.hh"
#include "Boundary_Source.hh"
#include "Boundary_Source_Parser.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Constructive_Solid_Geometry_Parser.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Material.hh"
#include "Material_Parser.hh"
#include "Solid_Geometry.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "Weight_Function.hh"
#include "XML_Node.hh"

using namespace std;

vector<shared_ptr<Weight_Function> > get_weight_functions(XML_Node input_node)
{
    // Get angular discretization
    Angular_Discretization_Parser angular_parser;
    shared_ptr<Angular_Discretization> angular
        = angular_parser.parse_from_xml(input_node.get_child("angular_discretization"));

    // Get energy discretization
    Energy_Discretization_Parser energy_parser;
    shared_ptr<Energy_Discretization> energy
        = energy_parser.parse_from_xml(input_node.get_child("energy_discretization"));

    // Get boundary sources
    Boundary_Source_Parser boundary_parser(angular,
                                           energy);
    vector<shared_ptr<Boundary_Source> > boundary_sources
        = boundary_parser.parse_from_xml(input_file.get_child("boundary_sources"));

    // Get materials
    Material_Parser material_parser(angular,
                                    energy);
    vector<shared_ptr<Material> > materials
        = material_parser.parse_from_xml(input_file.get_child("materials"));

    // Get solid geometry
    Constructive_Solid_Geometry_Parser solid_parser(materials,
                                                    boundary_sources);
    shared_ptr<Constructive_Solid_Geometry> solid_geometry
        = solid_parser.parse_from_xml(input_node);
    int dimension = solid_geometry->dimension();

    // Parser for basis and weight functions
    Weak_Spatial_Discretization_Parser spatial_parser(solid_geometry);

    // Get basis functions
    XML_Node bases_node = input_node.get_child("basis_functions");
    int number_of_bases = bases_node.get_child("number_of_basis_functions");
    vector<shared_ptr<Basis_Function> > basis_functions
        = spatial_parser.get_basis_functions(input_node.get_child("basis_functions"),
                                             number_of_bases,
                                             dimension);

    // Get weight functions
    XML_Node weights_node = input_node.get_child("weight_functions");
    int number_of_weights = weights_node.get_child("number_of_weight_functions");
    vector<shared_ptr<Weight_Function> > weight_functions
        = spatial_parser.get_weight_functions(input_node.get_child("weight_functions"),
                                              number_of_weights,
                                              dimension,
                                              basis_functions);
}

int test_integrals(string input_filename)
{
    int checksum = 0;
    
    XML_Document input_file(input_filename);
    
    vector<shared_ptr<Weight_Function> > weight_functions
        = get_weight_functions(input_file.get_child("input"));

    return checksum;
}

int main(int argc, char **argv)
{
    int checksum = 0;

    if (argc != 2)
    {
        cerr << "usage: tst_Weight_Function [input_folder]" << endl;
        return 1;
    }
    
    string input_folder = argv[1];
    input_folder += "/";
    {
        string input_filename = input_folder + "/weight_function.xml";
        
        checksum += test_integrals(input_filename);
    }
    
    return checksum;
}
