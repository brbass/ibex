#include <functional>

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source_Parser.hh"
#include "Cartesian_Plane.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Interpolation_Solid_Geometry.hh"
#include "Material.hh"
#include "Material_Parser.hh"
#include "Solid_Geometry.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "Weight_Function.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;

shared_ptr<Weak_Spatial_Discretization> get_spatial(int dimension,
                                                    function<double(vector<double>)> const &source,
                                                    XML_Node input_node)
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

    // Get solid geometry
    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces(2 * dimension);
    for (int d = 0; d < dimension; ++d)
    {
        for (int i = 0; i < 2; ++i)
        {
            int index = d + dimension * i;
            double position = i == 0 ? -1 : 1;
            double normal = i == 0 ? -1 : 1;
            boundary_surfaces[index]
                = make_shared<Cartesian_Plane>(index,
                                               dimension,
                                               Surface::Surface_Type::BOUNDARY,
                                               d,
                                               position,
                                               normal);
            boundary_surfaces[index]->set_boundary_source(boundary_sources[0]);
        }
                                               
    }
    shared_ptr<Solid_Geometry> solid_geometry
        = make_shared<Analytic_Solid_Geometry>(dimension,
                                               angular,
                                               energy,
                                               source);
    
    // Parser for basis and weight functions
    Weak_Spatial_Discretization_Parser spatial_parser(solid_geometry,
                                                      boundary_surfaces);
    return spatial_parser.get_weak_discretization(input_node.get_child("spatial_discretization"));
}

int test_interpolation(int dimension,
                       function<double(vector<double>)> const &source,
                       function<double(vector<double>)> const &d_source,
                       XML_Node input_node)
{
    int checksum = 0;
    
    shared_ptr<Weak_Spatial_Discretization> spatial
        = get_spatial(dimension,
                      source,
                      input_node);
    
    return checksum;
}

int main(int argc, char **argv)
{
    int checksum = 0;

    if (argc != 2)
    {
        cerr << "usage: tst_Interpolation [input_folder]" << endl;
        return 1;
    }
    
    string input_folder = argv[1];
    input_folder += "/";
    vector<string> input_filenames
        = {input_folder + "mls_interpolation.xml"};
    
    for (string input_filename : input_filenames)
    {
        // Get XML document
        XML_Document input_file(input_filename);
        XML_Node input_node = input_file.get_child("input");
        int dimension = input_node.get_child("angular_discretization").get_child_value<int>("dimension");
        
        // Test constant
        {
            function<double(vector<double>)> source
                = [](vector<double> const &position)
                {
                    return 1.;
                };
            function<double(vector<double>)> d_source
                = [](vector<double> const &position)
                {
                    return 0.;
                };
            
            checksum += test_interpolation(dimension,
                                           source,
                                           d_source,
                                           input_node);
        }
    }
    
    return checksum;
}
