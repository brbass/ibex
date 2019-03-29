#include <cmath>
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <mpi.h>
#if defined(ENABLE_OPENMP)
    #include <omp.h>
#else
    inline void omp_set_num_threads(int i) {return;}
#endif

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source.hh"
#include "Boundary_Source_Parser.hh"
#include "Cartesian_Plane.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Constructive_Solid_Geometry_Parser.hh"
#include "Cross_Section.hh"
#include "Discrete_To_Moment.hh"
#include "Discrete_Value_Operator.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Material.hh"
#include "Material_Parser.hh"
#include "Meshless_Sweep.hh"
#include "Meshless_Sweep_Parser.hh"
#include "Moment_Value_Operator.hh"
#include "SUPG_Internal_Source_Operator.hh"
#include "SUPG_Moment_To_Discrete.hh"
#include "Transport_Discretization.hh"
#include "Vector_Functions.hh"
#include "Vector_Operator_Functions.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;



double get_solution(std::shared_ptr<Constructive_Solid_Geometry> solid,
                    std::shared_ptr<Angular_Discretization> angular,
                    std::vector<double> const &position,
                    std::vector<double> const &direction,
                    bool &edge,
                    bool &corner)
{
    // Reverse direction
    std::vector<double> reverse_direction = {-direction[0], -direction[1]};

    // Get boundary surface index and distance
    vector<int> all_boundary_surfaces = solid->find_all_boundary_surfaces(position);
    int const region = 0;
    int boundary;
    double distance;

    edge = all_boundary_surfaces.size() > 0 ? true : false;
    corner = all_boundary_surfaces.size() > 1 ? true : false;
    
    if (all_boundary_surfaces.size() == 0) // internal point
    {
        int boundary_region;
        vector<double> final_position;
        boundary = solid->next_boundary(region, // initial region
                                        position,
                                        reverse_direction,
                                        boundary_region,
                                        distance,
                                        final_position);
    }
    else // point on at least one boundary
    {
        // cout << position[0] << "\t" << position[1] << "\t" << direction[0] << "\t" << direction[1];
        // If point is inward normal for any boundary, assign boundary and distance
        bool inward = false;
        for (int surface_index : all_boundary_surfaces)
        {
            Surface::Normal normal = solid->cartesian_boundary_surface(surface_index)->normal_direction(position,
                                                                                                        false); // check point
            if (Vector_Functions::dot(normal.direction, direction) < 0)
            {
                // cout << "\tin\t";
                // cout << normal.direction[0] << "\t" << normal.direction[1];
                boundary = surface_index;
                distance = 0;
                inward = true;
            }
        }
        
        // If point is outward normal for every point, can calculate analytic solution
        if (!inward)
        {
            // cout << "\tall out";
            // Go a small distance backwards
            double delta = solid->delta_distance();
            vector<double> inside_position;
            solid->new_position(delta,
                                position,
                                reverse_direction,
                                inside_position);

            // Get point information
            int boundary_region;
            vector<double> final_position;
            boundary = solid->next_boundary(region,
                                            inside_position,
                                            reverse_direction,
                                            boundary_region,
                                            distance,
                                            final_position);
            distance += delta;
        }
        // cout << endl;
    }
    Assert(boundary != Solid_Geometry::Geometry_Errors::NO_SURFACE);

    // Get boundary surface
    std::shared_ptr<Cartesian_Plane> boundary_surface = solid->cartesian_boundary_surface(boundary);
    
    // Get material and boundary source
    int material_index = 0;
    std::shared_ptr<Material> material = solid->material(material_index);
    std::shared_ptr<Boundary_Source> boundary_source = boundary_surface->boundary_source();
    Assert(material);
    Assert(boundary_source);
    
    // Get physical data : assumes isotropic boundary source
    double const psi0 = boundary_source->data()[0];
    double const sigma_t = material->sigma_t()->data()[0];
    double const q = material->internal_source()->data()[0];

    // Solve for expected angular flux value
    double const k = exp(-sigma_t * distance);
    double const d = angular->angular_normalization();
    return psi0 * k + q / (sigma_t * d) * (1 - k);
}

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
    for (shared_ptr<Boundary_Source> source : boundary_sources)
    {
        Assert(source->dependencies().angular == Boundary_Source::Dependencies::Angular::ISOTROPIC);
    }
    
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

    // Get discrete to moment operator
    shared_ptr<Discrete_To_Moment> D
        = make_shared<Discrete_To_Moment>(spatial,
                                          angular,
                                          energy);
    
    // Get discrete and moment value operators
    shared_ptr<Discrete_Value_Operator> discrete_value
        = make_shared<Discrete_Value_Operator>(spatial,
                                               angular,
                                               energy,
                                               false); // weighted
    shared_ptr<Moment_Value_Operator> moment_value
        = make_shared<Moment_Value_Operator>(spatial,
                                             angular,
                                             energy,
                                             false); // weighted
    
    // Get combined operator
    shared_ptr<Vector_Operator> combined = Linv * M * Q;
    
    // Get solution vector
    vector<double> coefficients(Q->column_size());

    // Get solution coefficients
    (*combined)(coefficients);

    // Get moment coefficients
    vector<double> moment_coefficients = coefficients;
    (*D)(moment_coefficients);

    // Get solution
    vector<double> solution = coefficients;
    (*discrete_value)(solution);
    vector<double> moment_solution = moment_coefficients;
    (*moment_value)(moment_solution);

    // Check solution
    int number_of_points = spatial->number_of_points();
    int number_of_ordinates = angular->number_of_ordinates();
    int number_of_moments = angular->number_of_moments();
    double err_phi = 0.;
    double err_psi = 0.;
    double sum_phi = 0.;
    double sum_psi = 0.;
    vector<double> const weights = angular->weights();
    for (int i = 0; i < number_of_points; ++i)
    {
        // Get position
        double sol_num = 0.;
        double sol_ana = 0.;
        vector<double> const position = spatial->point(i)->position();
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            // Get direction
            std::vector<double> const direction = angular->direction(o);
            
            // Get analytic solution
            bool edge;
            bool corner;
            double const analytic = get_solution(solid,
                                                 angular,
                                                 position,
                                                 direction,
                                                 edge,
                                                 corner);

            // Compare numeric to analytic solution
            int k = o + number_of_ordinates * i;

            double const weight = weights[o];
            sol_ana += analytic * weight;
            sum_psi += analytic;
            err_psi += abs(analytic - solution[k]);
            // cout << "ana: " << analytic << "\t" << "num: " << solution[k] << "\t" << "err: " << analytic - solution[k] << endl;
        }
        int k = number_of_moments * i;
        sum_phi += sol_ana;
        err_phi += abs(moment_solution[k] - sol_ana);
        // for (int d = 0; d < solid->dimension(); ++d)
        // {
        //     cout << position[d] << "\t";
        // }
        // cout << "ana: " << sol_ana << "\t" << "num: " << moment_solution[k] << "\t" << "err: " << sol_ana - moment_solution[k] << endl;
    }
    err_phi /= sum_phi;
    err_psi /= sum_psi;
    cout << "err_psi: " << err_psi << "\t" << "err_phi: " << err_phi << endl;
}

int main(int argc, char **argv)
{
    // MPI init
    MPI_Init(&argc, &argv);
    
    // Get input file
    string input_filename = argv[1];
    XML_Document input_file(input_filename);
    XML_Node input_node = input_file.get_child("input");

    // Set number of procs
    int number_of_threads = input_node.get_attribute<int>("number_of_threads");
    omp_set_num_threads(number_of_threads);

    // Get output file
    string output_filename = input_filename + ".out";
    XML_Document output_file;
    XML_Node output_node = output_file.append_child("output");

    // Run the problem
    run_problem(input_node,
                output_node);
    output_file.save(output_filename);

    // MPI finalize
    MPI_Finalize();
}
