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
                   shared_ptr<Constructive_Solid_Geometry> &solid,
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
    solid
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

double get_solution(double sigma_t,
                    double angular_normalization,
                    double internal_source,
                    double boundary_source,
                    double distance)
{
    double att = exp(-sigma_t * distance);

    return boundary_source * att + internal_source / (angular_normalization * sigma_t) * (1 - att);
}

int run_test(string input_filename)
{
    // Parse input file
    shared_ptr<Weak_Spatial_Discretization> spatial;
    shared_ptr<Angular_Discretization> angular;
    shared_ptr<Energy_Discretization> energy;
    shared_ptr<Transport_Discretization> transport;
    shared_ptr<Constructive_Solid_Geometry> solid;
    vector<shared_ptr<Material> > materials;
    vector<shared_ptr<Boundary_Source> > sources;
    shared_ptr<Weak_RBF_Sweep> sweeper;
    get_transport(input_filename,
                  spatial,
                  angular,
                  energy,
                  transport,
                  solid,
                  materials,
                  sources,
                  sweeper);
    shared_ptr<Material> material = materials[0];
    shared_ptr<Boundary_Source> source = sources[0];
    
    // Get data
    int phi_size = transport->phi_size();
    int psi_size = transport->psi_size();
    int dimension = solid->dimension();
    int number_of_augments = transport->number_of_augments();
    int number_of_points = spatial->number_of_points();
    int number_of_groups = energy->number_of_groups();
    int number_of_ordinates = angular->number_of_ordinates();
    double angular_normalization = angular->angular_normalization();
    
    // Get RHS
    vector<double> rhs(psi_size + number_of_augments, 0);
    for (int i = 0; i < number_of_points; ++i)
    {
        vector<double> const internal_source = spatial->weight(i)->material()->internal_source()->data();
        for (int g = 0; g < number_of_groups; ++g)
        {
            double is = internal_source[g];
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                int k = g + number_of_groups * (o + number_of_ordinates * i);
                rhs[k] = is / angular_normalization;
            }
        }
    }
    
    // Sweep
    sweeper->set_include_boundary_source(true);
    (*sweeper)(rhs);
    
    // Get values
    shared_ptr<Discrete_Value_Operator> discrete_value
        = make_shared<Discrete_Value_Operator>(spatial,
                                               angular,
                                               energy,
                                               false); // weighted
    (*discrete_value)(rhs);
    
    // Check results
    vector<double> const internal_source = material->internal_source()->data();
    vector<double> const sigma_t = material->sigma_t()->data();
    vector<double> const boundary_source = source->data();
    for (int i = 0; i < number_of_points; ++i)
    {
        vector<double> const position = spatial->weight(i)->position();
        
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            // Get distance to boundary
            vector<double> const direction = angular->direction(o);
            vector<double> opp_direction(dimension);
            for (int d = 0; d < dimension; ++d)
            {
                opp_direction[d] = -direction[d];
            }
            int boundary_region;
            double distance;
            vector<double> final_position;
            solid->next_boundary(position,
                                 opp_direction,
                                 boundary_region,
                                 distance,
                                 final_position);
            
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k = g + number_of_groups * (o + number_of_ordinates * i);
                double solution = get_solution(sigma_t[g],
                                               angular_normalization,
                                               internal_source[g],
                                               boundary_source[g + number_of_groups * o],
                                               distance);
                
                cout << rhs[k] << "\t" << solution << endl;
            }
        }
    }

    return 0;
}

int main(int argc, char **argv)
{
    int checksum = 0;
    
    MPI_Init(&argc, &argv);
    
    if (argc != 2)
    {
        cerr << "usage: tst_Purely_Absorbing [input_folder]" << endl;
        return 1;
    }
    
    string input_folder = argv[1]; 
    input_folder += "/";
    {
        string input_filename = input_folder + "boundary_slab.xml";
        
        checksum += run_test(input_filename);
    }
    
    MPI_Finalize();
    
    return checksum;
}
