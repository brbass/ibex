#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Parser.hh"
#include "Basis_Function.hh"
#include "Boundary_Source.hh"
#include "Boundary_Source_Parser.hh"
#include "Check_Equality.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Constructive_Solid_Geometry_Parser.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Material.hh"
#include "Material_Parser.hh"
#include "Solid_Geometry.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "Weight_Function.hh"
#include "XML_Document.hh"

using namespace std;
namespace ce = Check_Equality;

int check_results(vector<double> const &vec,
                   vector<double> const &ana,
                   string description,
                   int index,
                   double tolerance)
{
    int checksum = 0;
    int w = 15;
    int size = vec.size();
    
    cout << description << " results for weight (" << index << ")" << endl;
    if (!ce::approx(vec, ana, tolerance))
    {
        cout << "FAILED" << endl;
        checksum += 1;
    }
    cout << setw(w) << "calculated";
    cout << setw(w) << "expected";
    cout << setw(w) << "error";
    cout << endl;
    for (int i = 0; i < size; ++i)
    {
        cout << setw(w) << vec[i];
        cout << setw(w) << ana[i];
        cout << setw(w) << vec[i] - ana[i];
        cout << endl;
    }
    cout << endl;
    
    return checksum;
}

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
        = boundary_parser.parse_from_xml(input_node.get_child("boundary_sources"));

    // Get materials
    Material_Parser material_parser(angular,
                                    energy);
    vector<shared_ptr<Material> > materials
        = material_parser.parse_from_xml(input_node.get_child("materials"));

    // Get solid geometry
    Constructive_Solid_Geometry_Parser solid_parser(materials,
                                                    boundary_sources);
    shared_ptr<Constructive_Solid_Geometry> solid_geometry
        = solid_parser.parse_from_xml(input_node.get_child("solid_geometry"));
    int dimension = solid_geometry->dimension();

    // Parser for basis and weight functions
    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces
        = solid_geometry->cartesian_boundary_surfaces();
    Weak_Spatial_Discretization_Parser spatial_parser(solid_geometry,
                                                      boundary_surfaces);
    
    // Get basis functions
    XML_Node bases_node = input_node.get_child("basis_functions");
    int number_of_bases = bases_node.get_child_value<int>("number_of_basis_functions");
    vector<shared_ptr<Basis_Function> > basis_functions
        = spatial_parser.get_basis_functions(bases_node,
                                             number_of_bases,
                                             dimension);

    // Get weight functions
    XML_Node weights_node = input_node.get_child("weight_functions");
    int number_of_weights = weights_node.get_child_value<int>("number_of_weight_functions");
    return spatial_parser.get_weight_functions(weights_node,
                                               number_of_weights,
                                               dimension,
                                               basis_functions);
}

int test_integrals(string input_filename)
{
    int checksum = 0;

    // Get XML document
    XML_Document input_file(input_filename);
    XML_Node input_node = input_file.get_child("input");

    // Get weight functions
    vector<shared_ptr<Weight_Function> > weight_functions
        = get_weight_functions(input_node);

    // Loop through weight functions to compare results
    XML_Node results_node = input_node.get_child("expected_integrals");
    for (XML_Node node = results_node.get_child("weight");
         node;
         node = node.get_sibling("weight",
                                 false))
    {
        int index = node.get_attribute<int>("index");

        // Get weight function information
        shared_ptr<Weight_Function> weight = weight_functions[index];
        Weight_Function::Integrals const integrals = weight->integrals();
        int number_of_basis_functions = weight->number_of_basis_functions();
        int dimension = weight->dimension();
        int number_of_boundary_surfaces = weight->number_of_boundary_surfaces();

        // Check volume results
        vector<double> const &iv_w = integrals.iv_w;
        vector<double> const ana_iv_w = node.get_child_vector<double>("iv_w", 1);
        checksum += check_results(iv_w,
                                  ana_iv_w,
                                  "iv_w",
                                  index,
                                  1e-7);
        
        vector<double> const &iv_dw = integrals.iv_dw;
        vector<double> const ana_iv_dw = node.get_child_vector<double>("iv_dw", dimension);
        checksum += check_results(iv_dw,
                                  ana_iv_dw,
                                  "iv_dw",
                                  index,
                                  1e-7);

        vector<double> const &iv_b_w = integrals.iv_b_w;
        vector<double> const ana_iv_b_w = node.get_child_vector<double>("iv_b_w", number_of_basis_functions);
        checksum += check_results(iv_b_w,
                                  ana_iv_b_w,
                                  "iv_b_w",
                                  index,
                                  1e-7);

        vector<double> const &iv_b_dw = integrals.iv_b_dw;
        vector<double> const ana_iv_b_dw = node.get_child_vector<double>("iv_b_dw", number_of_basis_functions * dimension);
        checksum += check_results(iv_b_dw,
                                  ana_iv_b_dw,
                                  "iv_b_dw",
                                  index,
                                  1e-7);

        vector<double> const &iv_db_w = integrals.iv_db_w;
        vector<double> const ana_iv_db_w = node.get_child_vector<double>("iv_db_w", number_of_basis_functions * dimension);
        checksum += check_results(iv_db_w,
                                  ana_iv_db_w,
                                  "iv_db_w",
                                  index,
                                  1e-7);

        vector<double> const &iv_db_dw = integrals.iv_db_dw;
        vector<double> const ana_iv_db_dw = node.get_child_vector<double>("iv_db_dw", number_of_basis_functions * dimension * dimension);
        checksum += check_results(iv_db_dw,
                                  ana_iv_db_dw,
                                  "iv_db_dw",
                                  index,
                                  1e-7);
        
        // Check surface results
        if (number_of_boundary_surfaces > 0)
        {
            vector<double> const &is_w = integrals.is_w;
            vector<double> const ana_is_w = node.get_child_vector<double>("is_w", number_of_boundary_surfaces);
            checksum += check_results(is_w,
                                      ana_is_w,
                                      "is_w",
                                      index,
                                      1e-12);
            
            vector<double> const &is_b_w = integrals.is_b_w;
            vector<double> const ana_is_b_w = node.get_child_vector<double>("is_b_w", number_of_boundary_surfaces * number_of_basis_functions);
            checksum += check_results(is_b_w,
                                      ana_is_b_w,
                                      "is_b_w",
                                      index,
                                      1e-12);
        }
    }
    
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
