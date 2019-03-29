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
#include "Arbitrary_Discrete_Value_Operator.hh"
#include "Arbitrary_Moment_Value_Operator.hh"
#include "Boundary_Source.hh"
#include "Boundary_Source_Parser.hh"
#include "Cartesian_Plane.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Constructive_Solid_Geometry_Parser.hh"
#include "Conversion.hh"
#include "Cross_Section.hh"
#include "Discrete_To_Moment.hh"
#include "Discrete_Value_Operator.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Integral_Error_Operator.hh"
#include "Integration_Mesh.hh"
#include "Material.hh"
#include "Material_Parser.hh"
#include "Meshless_Function_Factory.hh"
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

    // Get combined operator
    shared_ptr<Vector_Operator> combined = Linv * M * Q;
    
    // Get solution vector
    vector<double> coefficients(Q->column_size());

    // Get solution coefficients
    (*combined)(coefficients);

    // Get moment coefficients
    vector<double> moment_coefficients = coefficients;
    (*D)(moment_coefficients);

    // Get solution functions
    int number_of_ordinates = angular->number_of_ordinates();
    int number_of_moments = angular->number_of_moments();
    std::function<vector<double>(vector<double>)> moment_solution
        = [&](vector<double> const& position)
        {
            vector<double> result(number_of_ordinates);
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                std::vector<double> const direction = angular->direction(o);
                bool corner;
                bool edge;
                result[o] = get_solution(solid,
                                         angular,
                                         position,
                                         direction,
                                         edge,
                                         corner);
            }
            angular->discrete_to_moment(result);
            return result;
        };

    std::function<vector<double>(vector<double>)> discrete_solution
        = [&](vector<double> const& position)
        {
            vector<double> result(number_of_ordinates);
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                std::vector<double> const direction = angular->direction(o);
                bool corner;
                bool edge;
                result[o] = get_solution(solid,
                                         angular,
                                         position,
                                         direction,
                                         edge,
                                         corner);
            }
            return result;
        };
    
    // Get integral errors
    vector<string> error_types;
    vector<vector<double> > errors;
    XML_Node errors_node = input_node.get_child("errors");
    for (XML_Node error_node = errors_node.get_child("error");
         error_node;
         error_node = error_node.get_sibling("error",
                                             false))
    {
        // Get integration options
        int dimension = spatial->dimension();
        std::shared_ptr<Integration_Mesh_Options> integration_options
            = make_shared<Integration_Mesh_Options>();
        integration_options->initialize_from_weak_options(spatial->options());
        integration_options->adaptive_quadrature
            = error_node.get_attribute<bool>("adaptive_quadrature",
                                             false);
        if (integration_options->adaptive_quadrature)
        {
            integration_options->minimum_radius_ordinates
                = error_node.get_attribute<bool>("minimum_radius_ordinates");
        }
        
        integration_options->integration_ordinates
            = error_node.get_attribute<int>("integration_ordinates",
                                            integration_options->integration_ordinates);
        integration_options->limits
            = error_node.get_child_matrix<double>("limits",
                                                  dimension,
                                                  2,
                                                  integration_options->limits);
        integration_options->dimensional_cells
            = error_node.get_child_vector<int>("dimensional_cells",
                                               dimension,
                                               integration_options->dimensional_cells);

        // Get error operator options
        Integral_Error_Operator::Options options;
        string norm_str = error_node.get_attribute<string>("norm");
        string angular_str = error_node.get_attribute<string>("angular");
        string energy_str = "group";
        options.norm = options.norm_conversion()->convert(norm_str);
        options.angular = options.angular_conversion()->convert(angular_str);
        options.energy = options.energy_conversion()->convert(energy_str);

        // Get analytic solution
        std::function<vector<double>(vector<double>)> solution;
        std::vector<double> error;
        switch (options.angular)
        {
        case Integral_Error_Operator::Options::Angular::MOMENTS:
            solution = moment_solution;
            error = moment_coefficients;
            break;
        case Integral_Error_Operator::Options::Angular::ORDINATES:
            solution = discrete_solution;
            error = coefficients;
            break;
        default:
            AssertMsg(false, "integral error angular type " + angular_str + "not found");
        }
        
        shared_ptr<Integral_Error_Operator> op
            = make_shared<Integral_Error_Operator>(options,
                                                   integration_options,
                                                   spatial,
                                                   angular,
                                                   energy,
                                                   solution);

        (*op)(error);
        error_types.push_back(norm_str + "_" + angular_str);
        errors.push_back(error);
    }

    vector<string> value_types;
    vector<vector<double> > values;
    vector<vector<vector<double>>> value_points;
    XML_Node values_node = input_node.get_child("values");
    for (XML_Node value_node = values_node.get_child("value");
         value_node;
         value_node = value_node.get_sibling("value",
                                             false))
    {
        // Get spatial data
        int dimension = spatial->dimension();
        vector<vector<double> > limits = spatial->options()->limits;

        // Get number of points
        vector<int> dimensional_points
            = value_node.get_child_vector<int>("points",
                                               dimension);

        // Get Cartesian points
        Meshless_Function_Factory factory;
        int num_eval_points;
        vector<vector<double> > eval_points;
        factory.get_cartesian_points(dimension,
                                     dimensional_points,
                                     limits,
                                     num_eval_points,
                                     eval_points);

        // Evaluate values
        vector<double> value;
        string eval_type = value_node.get_attribute<string>("type");
        if (eval_type == "moments")
        {
            shared_ptr<Arbitrary_Moment_Value_Operator> op
                = make_shared<Arbitrary_Moment_Value_Operator>(spatial,
                                                               angular,
                                                               energy,
                                                               eval_points);
            value = moment_coefficients;
            (*op)(value);
        }
        else if (eval_type == "ordinates")
        {
            shared_ptr<Arbitrary_Discrete_Value_Operator> op
                = make_shared<Arbitrary_Discrete_Value_Operator>(spatial,
                                                                 angular,
                                                                 energy,
                                                                 eval_points);
            value = coefficients;
            (*op)(value); 
        }
        else if (eval_type == "analytic_moments")
        {
            value.resize(num_eval_points);
            for (int i = 0; i < num_eval_points; ++i)
            {
                value[i] = moment_solution(eval_points[i])[0];
            }
        }
        else if (eval_type == "analytic_ordinates")
        {
            value.resize(num_eval_points * number_of_ordinates);
            for (int i = 0; i < num_eval_points; ++i)
            {
                vector<double> sol = discrete_solution(eval_points[i]);
                for (int o = 0; o < number_of_ordinates; ++o)
                {
                    value[o + number_of_ordinates * i] = sol[o];
                }
            }
        }
        else
        {
            cout << "evaluation type " << eval_type << " not found and thus ignored" << endl;
        }
        values.push_back(value);
        value_types.push_back(eval_type);
        value_points.push_back(eval_points);
    }
    
    energy->output(output_node.append_child("energy_discretization"));
    angular->output(output_node.append_child("angular_discretization"));
    spatial->output(output_node.append_child("spatial_discretization"));
    transport->output(output_node.append_child("transport_discretization"));
    solid->output(output_node.append_child("solid_geometry"));
    Linv->output(output_node.append_child("transport"));
    if (errors.size() > 0)
    {
        XML_Node error_node = output_node.append_child("error");
        for (int i = 0; i < errors.size(); ++i)
        {
            error_node.set_child_vector<double>(errors[i], error_types[i], "angular-cell");
        }
    }
    if (values.size() > 0)
    {
        XML_Node values_node = output_node.append_child("values");
        for (int i = 0; i < values.size(); ++i)
        {
            XML_Node value_node = values_node.append_child("value");
            value_node.set_attribute(value_types[i], "type");
            value_node.set_child_matrix<double>(value_points[i], "points");
            value_node.set_child_vector<double>(values[i], "result", "angular-cell");
        }
    }
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
