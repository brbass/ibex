#include <mpi.h>

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Factory.hh"
#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source.hh"
#include "Boundary_Source_Parser.hh"
#include "Cartesian_Plane.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Constructive_Solid_Geometry_Parser.hh"
#include "Cross_Section.hh"
#include "Dimensional_Moments.hh"
#include "Discrete_Value_Operator.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Material.hh"
#include "Material_Factory.hh"
#include "Material_Parser.hh"
#include "Region.hh"
#include "Transport_Discretization.hh"
#include "Weak_Meshless_Sweep.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Factory.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

#include <iomanip>
#include <iostream>
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

void get_one_region(bool basis_mls,
                    bool weight_mls,
                    string basis_type,
                    string weight_type,
                    shared_ptr<Weight_Function_Options> weight_options,
                    shared_ptr<Weak_Spatial_Discretization_Options> weak_options,
                    int dimension,
                    int angular_rule,
                    int num_dimensional_points,
                    double radius_num_intervals,
                    double sigma_t,
                    double internal_source,
                    double boundary_source,
                    double length,
                    shared_ptr<Weak_Spatial_Discretization> &spatial,
                    shared_ptr<Angular_Discretization> &angular,
                    shared_ptr<Energy_Discretization> &energy,
                    shared_ptr<Transport_Discretization> &transport,
                    shared_ptr<Constructive_Solid_Geometry> &solid,
                    vector<shared_ptr<Material> > &materials,
                    vector<shared_ptr<Boundary_Source> > &boundary_sources,
                    shared_ptr<Meshless_Sweep> &sweeper)
{
    // Get angular discretization
    int number_of_moments = 1;
    Angular_Discretization_Factory angular_factory;
    angular = angular_factory.get_angular_discretization(dimension,
                                                         number_of_moments,
                                                         angular_rule);
    
    // Get energy discretization
    int number_of_groups = 1;
    energy = make_shared<Energy_Discretization>(number_of_groups);
    
    // Get material
    materials.resize(1);
    Material_Factory material_factory(angular,
                                      energy);
    materials[0]
        = material_factory.get_standard_material(0, // index
                                                 {sigma_t},
                                                 {0}, // sigma_s
                                                 {0}, // nu
                                                 {0}, // sigma_f
                                                 {0}, // chi
                                                 {internal_source});
    
    // Get boundary source
    boundary_sources.resize(1);
    Boundary_Source::Dependencies boundary_dependencies;
    boundary_sources[0]
        = make_shared<Boundary_Source>(0, // index
                                       boundary_dependencies,
                                       angular,
                                       energy,
                                       vector<double>(1, boundary_source),
                                       vector<double>(1, 0));
    
    // Get solid geometry
    vector<shared_ptr<Surface> > surfaces(2 * dimension);
    vector<shared_ptr<Region> > regions(1);
    for (int d = 0; d < dimension; ++d)
    {
        int index1 = 0 + 2 * d;
        surfaces[index1]
            = make_shared<Cartesian_Plane>(index1,
                                           dimension,
                                           Surface::Surface_Type::BOUNDARY,
                                           d,
                                           -0.5 * length,
                                           -1);
        int index2 = 1 + 2 * d;
        surfaces[index2]
            = make_shared<Cartesian_Plane>(index2,
                                           dimension,
                                           Surface::Surface_Type::BOUNDARY,
                                           d,
                                           0.5 * length,
                                           1);
    }
    for (shared_ptr<Surface> surface : surfaces)
    {
        surface->set_boundary_source(boundary_sources[0]);
    }
    vector<Surface::Relation> surface_relations(2 * dimension,
                                                Surface::Relation::NEGATIVE);
    regions[0]
        = make_shared<Region>(0, // index
                              materials[0],
                              surface_relations,
                              surfaces);
    solid
        = make_shared<Constructive_Solid_Geometry>(dimension,
                                                   surfaces,
                                                   regions,
                                                   materials,
                                                   boundary_sources);
    
    // Get spatial discretization
    Weak_Spatial_Discretization_Factory spatial_factory(solid,
                                                        solid->cartesian_boundary_surfaces());
    spatial
        = spatial_factory.get_simple_discretization(num_dimensional_points,
                                                    radius_num_intervals,
                                                    basis_mls,
                                                    weight_mls,
                                                    basis_type,
                                                    weight_type,
                                                    weight_options,
                                                    weak_options);
    
    // Get transport discretization
    transport
        = make_shared<Transport_Discretization>(spatial,
                                                angular,
                                                energy);
    
    // Get weak RBF sweep
    Meshless_Sweep::Options options;
    sweeper
        = make_shared<Weak_Meshless_Sweep>(options,
                                           spatial,
                                           angular,
                                           energy,
                                           transport);
}


void get_transport_from_xml(string input_filename,
                            shared_ptr<Weak_Spatial_Discretization> &spatial,
                            shared_ptr<Angular_Discretization> &angular,
                            shared_ptr<Energy_Discretization> &energy,
                            shared_ptr<Transport_Discretization> &transport,
                            shared_ptr<Constructive_Solid_Geometry> &solid,
                            vector<shared_ptr<Material> > &materials,
                            vector<shared_ptr<Boundary_Source> > &boundary_sources,
                            shared_ptr<Meshless_Sweep> &sweeper)
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
    Meshless_Sweep::Options options;
    sweeper
        = make_shared<Weak_Meshless_Sweep>(options,
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

int run_test(shared_ptr<Weak_Spatial_Discretization> spatial,
             shared_ptr<Angular_Discretization> angular,
             shared_ptr<Energy_Discretization> energy,
             shared_ptr<Transport_Discretization> transport,
             shared_ptr<Constructive_Solid_Geometry> solid,
             vector<shared_ptr<Material> > materials,
             vector<shared_ptr<Boundary_Source> > sources,
             shared_ptr<Meshless_Sweep> sweeper,
             bool print)
{
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
    int number_of_dimensional_moments = spatial->dimensional_moments()->number_of_dimensional_moments();
    double angular_normalization = 1;//angular->angular_normalization();
    
    // Get RHS
    vector<double> rhs(psi_size + number_of_augments, 0);
    for (int i = 0; i < number_of_points; ++i)
    {
        shared_ptr<Weight_Function> weight = spatial->weight(i);
        vector<double> const internal_source = weight->material()->internal_source()->data();
        double tau = weight->options()->tau;
        for (int g = 0; g < number_of_groups; ++g)
        {
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                vector<double> const direction = angular->direction(o);
                
                int k = g + number_of_groups * (o + number_of_ordinates * i);
                {
                    int d = 0;
                    rhs[k] = internal_source[d + number_of_dimensional_moments * g];
                }
                for (int d = 1; d < number_of_dimensional_moments; ++d)
                {
                    rhs[k] += tau * direction[d-1] * internal_source[d + number_of_dimensional_moments * g];
                }
                
                rhs[k] /= angular_normalization;
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
    int w = 16;
    int num_average = 0;
    double l2err = 0;
    double l2relerr = 0;
    if (print)
    {
        cout << setw(w) << "calculated" << setw(w) << "solution" << setw(w) << "rel_error" << endl;
    }
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
            int region;
            double distance;
            vector<double> final_position;
            solid->next_intersection(position,
                                     opp_direction,
                                     region,
                                     distance,
                                     final_position);

            // Next_intersection struggles on corners
            // Only include points that correctly find vacuum
            if (region == Solid_Geometry::Geometry_Errors::NO_REGION)
            {
                for (int g = 0; g < number_of_groups; ++g)
                {
                    int k = g + number_of_groups * (o + number_of_ordinates * i);
                
                    double solution = get_solution(sigma_t[g],
                                                   angular_normalization,
                                                   internal_source[g],
                                                   boundary_source[g + number_of_groups * o],
                                                   distance);
                    double err = rhs[k] - solution;
                    double relerr = err / solution;
                    l2err += err * err;
                    l2relerr += relerr * relerr;
                    num_average += 1;
                    if (o == 1 && print)
                    {
                        cout << setw(w) << rhs[k] << setw(w) << solution << setw(w) << relerr << endl;
                    }
                }
            }
        }
    }
    l2err = sqrt(l2err) / num_average;
    l2relerr = sqrt(l2relerr) / num_average;
    cout << setw(w) << l2err << setw(w) << l2relerr << endl;
    
    return 0;
}

int main(int argc, char **argv)
{
    int checksum = 0;
    
    MPI_Init(&argc, &argv);
    
    shared_ptr<Weak_Spatial_Discretization> spatial;
    shared_ptr<Angular_Discretization> angular;
    shared_ptr<Energy_Discretization> energy;
    shared_ptr<Transport_Discretization> transport;
    shared_ptr<Constructive_Solid_Geometry> solid;
    vector<shared_ptr<Material> > materials;
    vector<shared_ptr<Boundary_Source> > sources;
    shared_ptr<Meshless_Sweep> sweeper;
    bool print = false;
    
    if (argc == 2)
    {
        string input_folder = argv[1]; 
        input_folder += "/";
        string input_filename = input_folder + "boundary_slab.xml";
        
        get_transport_from_xml(input_filename,
                               spatial,
                               angular,
                               energy,
                               transport,
                               solid,
                               materials,
                               sources,
                               sweeper);
        
        checksum += run_test(spatial,
                             angular,
                             energy,
                             transport,
                             solid,
                             materials,
                             sources,
                             sweeper,
                             print);

    }
    else if (argc == 4)
    {
        int dimension = atoi(argv[1]);
        int angular_rule = dimension == 1 ? 2 : 1;
        int num_dimensional_points = atoi(argv[2]);
        double radius_num_intervals = atof(argv[3]);
        double sigma_t = 1.0;
        double internal_source = 2.0;
        double boundary_source = 1.0;
        double length = 2;
        bool basis_mls = true;
        bool weight_mls = true;
        string basis_type = "compact_gaussian";
        string weight_type = "compact_gaussian";
        shared_ptr<Weight_Function_Options> weight_options
            = make_shared<Weight_Function_Options>();
        shared_ptr<Weak_Spatial_Discretization_Options> weak_options
            = make_shared<Weak_Spatial_Discretization_Options>();
        weak_options->integration_ordinates = 64;
        weight_options->tau_const = 0.5;
        weak_options->include_supg = true;
        get_one_region(basis_mls,
                       weight_mls,
                       basis_type,
                       weight_type,
                       weight_options,
                       weak_options,
                       dimension,
                       angular_rule,
                       num_dimensional_points,
                       radius_num_intervals,
                       sigma_t,
                       internal_source,
                       boundary_source,
                       length,
                       spatial,
                       angular,
                       energy,
                       transport,
                       solid,
                       materials,
                       sources,
                       sweeper);

        cout << setw(16) << dimension << setw(16) << num_dimensional_points << setw(16) << radius_num_intervals;
        checksum += run_test(spatial,
                             angular,
                             energy,
                             transport,
                             solid,
                             materials,
                             sources,
                             sweeper,
                             print);
    }
    else
    {
        cerr << "usage: tst_Purely_Absorbing [input_folder]" << endl;
        cerr << "usage: tst_Purely_Absorbing [dimension num_points num_intervals]" << endl;
        return 1;
    }
    
    
    MPI_Finalize();
    
    return checksum;
}

