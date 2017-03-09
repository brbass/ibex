#include <mpi.h>

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source.hh"
#include "Boundary_Source_Parser.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Constructive_Solid_Geometry_Parser.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Material.hh"
#include "Material_Parser.hh"
#include "Transport_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "Weak_RBF_Sweep.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

#include <memory>
#include <vector>

using namespace std;

XML_Document get_xml_document(string input_filename)
{
    // Get input file
    string input_folder = input_filename.substr(0, input_filename.find_last_of("/\\") + 1);
    XML_Document input_file(input_filename);
    XML_Node input_node = input_file.get_child("input");
    XML_Node spatial_node = input_node.get_child("spatial_discretization");
    
    // Append external basis and weight function information
    string external_filename = input_folder + spatial_node.get_attribute<string>("file");
    XML_Document external_document(external_filename);
    XML_Node external_node = external_document.get_child("spatial_discretization");
    spatial_node.prepend_node(external_node.get_child("number_of_points"));
    spatial_node.get_child("basis_functions").append_all(external_node.get_child("basis_functions"));
    spatial_node.get_child("weight_functions").append_all(external_node.get_child("weight_functions"));

    // Append solid geometry information
    input_node.append_node(external_node.get_child("solid_geometry"));
    
    return input_file;
}

void get_transport(string input_filename,
                   shared_ptr<Weak_Spatial_Discretization> &spatial,
                   shared_ptr<Angular_Discretization> &angular,
                   shared_ptr<Energy_Discretization> &energy,
                   shared_ptr<Transport_Discretization> &transport,
                   vector<shared_ptr<Material> > &materials,
                   vector<shared_ptr<Boundary_Source> > &boundary_sources,
                   shared_ptr<Weak_RBF_Sweep> &sweeper)
{
    // Get input document
    XML_Document input_file = get_xml_document(input_filename);
    XML_Node input_node = input_file.get_child("input");
    
    // Get angular discretization
    Angular_Discretization_Parser angular_parser;
    angular
        = angular_parser.parse_from_xml(input_node.get_child("angular_discretization"));

    // Get energy discretization
    Energy_Discretization_Parser energy_parser;
    energy
        = energy_parser.parse_from_xml(input_node.get_child("energy_discretization"));

    // Get materials
    Material_Parser material_parser(angular,
                                    energy);
    materials
        = material_parser.parse_from_xml(input_node.get_child("materials"));
    
    // Get boundary sources
    Boundary_Source_Parser boundary_parser(angular,
                                           energy);
    boundary_sources
        = boundary_parser.parse_from_xml(input_node.get_child("boundary_sources"));

    // Get solid geometry
    Constructive_Solid_Geometry_Parser solid_parser(materials,
                                                    boundary_sources);
    shared_ptr<Constructive_Solid_Geometry> solid
        = solid_parser.parse_from_xml(input_node.get_child("solid_geometry"));

    // Get spatial discretization
    Weak_Spatial_Discretization_Parser spatial_parser(solid,
                                                      solid->cartesian_boundary_surfaces());
    spatial
        = spatial_parser.get_weak_discretization(input_node.get_child("spatial_discretization"));

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
}

int test_boundary_slab(string input_folder)
{
    // Parse input file
    string input_filename = input_folder + "boundary_slab.xml";
    shared_ptr<Weak_Spatial_Discretization> spatial;
    shared_ptr<Angular_Discretization> angular;
    shared_ptr<Energy_Discretization> energy;
    shared_ptr<Transport_Discretization> transport;
    vector<shared_ptr<Material> > materials;
    vector<shared_ptr<Boundary_Source> > sources;
    shared_ptr<Weak_RBF_Sweep> sweeper;
    get_transport(input_filename,
                  spatial,
                  angular,
                  energy,
                  transport,
                  materials,
                  sources,
                  sweeper);
    shared_ptr<Material> material = materials[0];
    shared_ptr<Boundary_Source> source = sources[0];
    
    // Get RHS
    int psi_size = transport->psi_size();
    int number_of_augments = transport->number_of_augments();
    vector<double> rhs(psi_size + number_of_augments, 0);

    // Sweep
    sweeper->set_include_boundary_source(true);
    (*sweeper)(rhs);
    
    // Check results
    int number_of_points = spatial->number_of_points();
    int number_of_groups = energy->number_of_groups();
    int number_of_ordinates = angular->number_of_ordinates();
    for (int i = 0; i < number_of_points; ++i)
    {
        int o = 1;
        int g = 0;
        int k = g + number_of_groups * (o + number_of_ordinates * i);
        
        cout << rhs[k] << endl;
    }

    return 0;
}

int main(int argc, char **argv)
{
    int checksum = 0;

    MPI_Init(&argc, &argv);
    
    if (argc != 2)
    {
        cerr << "usage: tst_Slab_Transport [input_folder]" << endl;
        return 1;
    }
    
    string input_folder = argv[1];
    input_folder += "/";
    
    checksum += test_boundary_slab(input_folder);
    
    MPI_Finalize();
    
    return checksum;
}

