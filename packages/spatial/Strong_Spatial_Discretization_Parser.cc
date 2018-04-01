#include "Strong_Spatial_Discretization_Parser.hh"

#include <limits>

#include "Angular_Discretization.hh"
#include "Basis_Function.hh"
#include "Cartesian_Distance.hh"
#include "Cartesian_Plane.hh"
#include "Compact_Gaussian_RBF.hh"
#include "Conversion.hh"
#include "Dimensional_Moments.hh"
#include "Energy_Discretization.hh"
#include "Linear_MLS_Function.hh"
#include "KD_Tree.hh"
#include "Meshless_Function.hh"
#include "Meshless_Function_Factory.hh"
#include "Quadratic_MLS_Function.hh"
#include "RBF_Factory.hh"
#include "RBF_Function.hh"
#include "RBF_Parser.hh"
#include "Solid_Geometry.hh"
#include "Truncated_Gaussian_RBF.hh"
#include "Strong_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Factory.hh"
#include "Weight_Function.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;

Strong_Spatial_Discretization_Parser::
Strong_Spatial_Discretization_Parser(shared_ptr<Solid_Geometry> solid_geometry,
                                   vector<shared_ptr<Cartesian_Plane> > const &boundary_surfaces):
    solid_geometry_(solid_geometry),
    boundary_surfaces_(boundary_surfaces)
{
    weight_radii_ = 1000 * numeric_limits<double>::epsilon();
}

shared_ptr<Strong_Spatial_Discretization> Strong_Spatial_Discretization_Parser::
get_strong_discretization(XML_Node input_node) const
{
    // Initialize data
    string input_format = input_node.get_attribute<string>("input_format");

    // Get discretization
    if (input_format == "points")
    {
        return get_points_discretization(input_node);
    }
    // else if (input_format == "cartesian")
    // {
    //     return get_cartesian_discretization(input_node);
    // }
    else
    {
        AssertMsg(false, "input format (" + input_format + ") not found");
        return shared_ptr<Strong_Spatial_Discretization>();
    }
}

shared_ptr<Strong_Spatial_Discretization> Strong_Spatial_Discretization_Parser::
get_points_discretization(XML_Node input_node) const
{
    RBF_Parser rbf_parser;
    Meshless_Function_Factory meshless_factory;
    Weak_Spatial_Discretization_Factory weak_factory(solid_geometry_,
                                                     boundary_surfaces_);
    
    // Get options
    shared_ptr<Weight_Function_Options> weight_options
        = get_weight_options(input_node.get_child("options"));
    shared_ptr<Weak_Spatial_Discretization_Options> weak_options
        = get_weak_options(input_node.get_child("options"));
    weak_options->identical_basis_functions
        = Weak_Spatial_Discretization_Options::Identical_Basis_Functions::FALSE;
    
    // Initialize data
    int dimension = solid_geometry_->dimension();
    XML_Node weights_node = input_node.get_child("weight_functions");
    
    // Get dimensional moments
    shared_ptr<Dimensional_Moments> dimensional_moments
        = make_shared<Dimensional_Moments>(weak_options->include_supg,
                                           dimension);
    
    // Get points
    string points_filename = input_node.get_attribute<string>("points_file");
    XML_Document points_doc(points_filename);
    XML_Node points_node = points_doc.get_child("input").get_child("spatial_discretization");
    int number_of_points = points_node.get_child_value<int>("number_of_points");
    vector<double> points_input = points_node.get_child_vector<double>("points", number_of_points * dimension);
    vector<vector<double> > points(number_of_points, vector<double>(dimension));
    for (int i = 0; i < number_of_points; ++i)
    {
        for (int d = 0; d < dimension; ++d)
        {
            int k = d + dimension * i;
            points[i][d] = points_input[k];
        }
    }
    
    // Make KD tree
    shared_ptr<KD_Tree> kd_tree
        = make_shared<KD_Tree>(dimension,
                               number_of_points,
                               points);
    
    // Find radii
    XML_Node radius_node = weights_node.get_child("radius_calculation");
    string radii_method = radius_node.get_attribute<string>("method");
    int number_of_neighbors = radius_node.get_child_value<int>("number_of_neighbors");
    double radius_multiplier = radius_node.get_child_value<double>("radius_multiplier");
    vector<double> radii;
    if (radii_method == "nearest")
    {
        meshless_factory.get_radii_nearest(kd_tree,
                                           dimension,
                                           number_of_points,
                                           number_of_neighbors,
                                           radius_multiplier,
                                           radii);
    }
    else if (radii_method == "coverage")
    {
        meshless_factory.get_radii_coverage(kd_tree,
                                            dimension,
                                            number_of_points,
                                            number_of_neighbors,
                                            radius_multiplier,
                                            radii);
    }
    else
    {
        AssertMsg(false, "radii method (" + radii_method + ") not found");
    }
    
    // Get RBF
    shared_ptr<RBF> rbf = rbf_parser.parse_from_xml(weights_node.get_child("meshless_function"));
    shared_ptr<Distance> distance = make_shared<Cartesian_Distance>(dimension);
    
    // Get independent meshless functions for basis
    vector<shared_ptr<Meshless_Function> > meshless_functions;
    meshless_factory.get_rbf_functions(number_of_points,
                                       radii,
                                       points,
                                       rbf,
                                       distance,
                                       meshless_functions);

    // Get MLS for basis (if applicable)
    string meshless_type = weights_node.get_child("meshless_function").get_attribute<string>("type");
    if (meshless_type == "rbf")
    {
        // Do nothing
    }
    else
    {
        // Find neighbors for MLS
        vector<vector<int> > neighbors;
        vector<vector<double> > squared_distances;
        meshless_factory.get_neighbors(kd_tree,
                                       rbf->range() == RBF::Range::GLOBAL,
                                       dimension,
                                       number_of_points,
                                       radii,
                                       radii,
                                       points,
                                       neighbors,
                                       squared_distances);
        Assert(meshless_factory.check_point_conditioning(number_of_points,
                                                         radii,
                                                         neighbors,
                                                         squared_distances));
        if (meshless_type == "linear_mls")
        {
            vector<shared_ptr<Meshless_Function> > mls_functions;
            meshless_factory.get_mls_functions(1,
                                               number_of_points,
                                               meshless_functions,
                                               neighbors,
                                               mls_functions);
            meshless_functions.swap(mls_functions);
        }
        else if (meshless_type == "quadratic_mls")
        {
            vector<shared_ptr<Meshless_Function> > mls_functions;
            meshless_factory.get_mls_functions(2,
                                               number_of_points,
                                               meshless_functions,
                                               neighbors,
                                               mls_functions);
            meshless_functions.swap(mls_functions);
        }
        else
        {
            AssertMsg(false, "meshless type (" + meshless_type + ") not found");
        }
    }

    // Get basis functions
    vector<shared_ptr<Basis_Function> > basis_functions;
    weak_factory.get_basis_functions(number_of_points,
                                     meshless_functions,
                                     basis_functions);

    // Get meshless functions for weight: use same distance as basis
    RBF_Factory rbf_factory;
    shared_ptr<RBF> weight_rbf = rbf_factory.get_rbf("wendland11");
    vector<shared_ptr<Meshless_Function> > weight_meshless_functions;
    vector<double> weight_radii(number_of_points, weight_radii_);
    meshless_factory.get_rbf_functions(number_of_points,
                                       weight_radii,
                                       points,
                                       weight_rbf,
                                       distance,
                                       weight_meshless_functions);

    // Get neighbors for weight functions
    vector<vector<int> > neighbors(number_of_points);
    vector<vector<double> > squared_distances(number_of_points);
    meshless_factory.get_point_neighbors(kd_tree,
                                         dimension,
                                         number_of_points,
                                         radii,
                                         points,
                                         neighbors,
                                         squared_distances);
    
    // Get weight functions
    vector<shared_ptr<Weight_Function> > weight_functions;
    weak_factory.get_weight_functions(number_of_points,
                                      weight_options,
                                      weak_options,
                                      dimensional_moments,
                                      neighbors,
                                      weight_meshless_functions,
                                      basis_functions,
                                      weight_functions);
    
    // Create spatial discretization
    return make_shared<Strong_Spatial_Discretization>(basis_functions,
                                                    weight_functions,
                                                    dimensional_moments,
                                                    weak_options);
}

// shared_ptr<Weak_Spatial_Discretization> Strong_Spatial_Discretization_Parser::
// get_cartesian_discretization(XML_Node input_node) const
// {
//     RBF_Parser rbf_parser;
//     Meshless_Function_Factory meshless_factory;
//     Weak_Spatial_Discretization_Factory weak_factory(solid_geometry_,
//                                                      boundary_surfaces_);
    
//     // Get options
//     shared_ptr<Weight_Function_Options> weight_options
//         = get_weight_options(input_node.get_child("options"));
//     shared_ptr<Weak_Spatial_Discretization_Options> weak_options
//         = get_weak_options(input_node.get_child("options"));
//     Assert(weak_options->identical_basis_functions
//            != Weak_Spatial_Discretization_Options::Identical_Basis_Functions::FALSE);
//     weak_options->identical_basis_functions
//         = Weak_Spatial_Discretization_Options::Identical_Basis_Functions::TRUE;
    
//     // Initialize data
//     int dimension = solid_geometry_->dimension();
//     XML_Node weights_node = input_node.get_child("weight_functions");
    
//     // Get dimensional moments
//     shared_ptr<Dimensional_Moments> dimensional_moments
//         = make_shared<Dimensional_Moments>(false, // do not include supg
//                                            dimension);
    
//     // Get points
//     vector<int> dimensional_points =
//         input_node.get_child_vector<int>("dimensional_points",
//                                          dimension);
//     int number_of_points;
//     vector<vector<double> > points;
//     meshless_factory.get_cartesian_points(dimension,
//                                           dimensional_points,
//                                           weak_options->limits,
//                                           number_of_points,
//                                           points);
    
//     // Make KD tree
//     shared_ptr<KD_Tree> kd_tree
//         = make_shared<KD_Tree>(dimension,
//                                number_of_points,
//                                points);
    
//     // Find radii
//     XML_Node radius_node = weights_node.get_child("radius_calculation");
//     string radii_method = radius_node.get_attribute<string>("method");
//     int number_of_neighbors = radius_node.get_child_value<int>("number_of_neighbors");
//     double radius_multiplier = radius_node.get_child_value<double>("radius_multiplier");
//     vector<double> radii;
//     if (radii_method == "nearest")
//     {
//         meshless_factory.get_radii_nearest(kd_tree,
//                                            dimension,
//                                            number_of_points,
//                                            number_of_neighbors,
//                                            radius_multiplier,
//                                            radii);
//     }
//     else if (radii_method == "coverage")
//     {
//         meshless_factory.get_radii_coverage(kd_tree,
//                                             dimension,
//                                             number_of_points,
//                                             number_of_neighbors,
//                                             radius_multiplier,
//                                             radii);
//     }
//     else
//     {
//         AssertMsg(false, "radii method (" + radii_method + ") not found");
//     }
    
//     // Find neighbors
//     vector<vector<int> > neighbors;
//     vector<vector<double> > squared_distances;
//     meshless_factory.get_neighbors(kd_tree,
//                                    dimension,
//                                    number_of_points,
//                                    radii,
//                                    radii,
//                                    points,
//                                    neighbors,
//                                    squared_distances);
//     Assert(meshless_factory.check_point_conditioning(number_of_points,
//                                                      radii,
//                                                      neighbors,
//                                                      squared_distances));
    
//     // Get RBF
//     shared_ptr<RBF> rbf = rbf_parser.parse_from_xml(weights_node.get_child("meshless_function"));
//     shared_ptr<Distance> distance = make_shared<Cartesian_Distance>(dimension);
    
//     // Get independent meshless functions
//     vector<shared_ptr<Meshless_Function> > meshless_functions;
//     meshless_factory.get_rbf_functions(number_of_points,
//                                        radii,
//                                        points,
//                                        rbf,
//                                        distance,
//                                        meshless_functions);

//     // Get MLS
//     string meshless_type = weights_node.get_child("meshless_function").get_attribute<string>("type");
//     if (meshless_type == "rbf")
//     {
//         // Do nothing
//     }
//     else if (meshless_type == "linear_mls")
//     {
//         vector<shared_ptr<Meshless_Function> > mls_functions;
//         meshless_factory.get_mls_functions(1,
//                                            number_of_points,
//                                            meshless_functions,
//                                            neighbors,
//                                            mls_functions);
//         meshless_functions.swap(mls_functions);
//     }
//     else if (meshless_type == "quadratic_mls")
//     {
//         vector<shared_ptr<Meshless_Function> > mls_functions;
//         meshless_factory.get_mls_functions(2,
//                                            number_of_points,
//                                            meshless_functions,
//                                            neighbors,
//                                            mls_functions);
//         meshless_functions.swap(mls_functions);
//     }
//     else
//     {
//         AssertMsg(false, "meshless type (" + meshless_type + ") not found");
//     }
    
//     // Get basis functions
//     vector<shared_ptr<Basis_Function> > basis_functions;
//     weak_factory.get_basis_functions(number_of_points,
//                                      meshless_functions,
//                                      basis_functions);

//     // Get weight functions
//     vector<shared_ptr<Weight_Function> > weight_functions;
//     weak_factory.get_weight_functions(number_of_points,
//                                       weight_options,
//                                       weak_options,
//                                       dimensional_moments,
//                                       neighbors,
//                                       meshless_functions,
//                                       basis_functions,
//                                       weight_functions);
    
//     // Create spatial discretization
//     return make_shared<Weak_Spatial_Discretization>(basis_functions,
//                                                     weight_functions,
//                                                     dimensional_moments,
//                                                     weak_options);
// }

shared_ptr<Weight_Function_Options> Strong_Spatial_Discretization_Parser::
get_weight_options(XML_Node input_node) const
{
    shared_ptr<Weight_Function_Options> options
        = make_shared<Weight_Function_Options>();
    
    options->tau_const = 0;
    options->output_material
        = input_node.get_attribute<bool>("output_material",
                                         options->output_material);
    options->output_integrals
        = input_node.get_attribute<bool>("output_integrals",
                                         options->output_integrals);
    
    return options;
}

shared_ptr<Weak_Spatial_Discretization_Options> Strong_Spatial_Discretization_Parser::
get_weak_options(XML_Node input_node) const
{
    int dimension = solid_geometry_->dimension();
    shared_ptr<Weak_Spatial_Discretization_Options> options
        = make_shared<Weak_Spatial_Discretization_Options>();
    options->discretization = Weak_Spatial_Discretization_Options::Discretization::STRONG;
    Meshless_Function_Factory meshless_factory;
    
    // Get weighting method
    string weighting_string = input_node.get_attribute<string>("weighting");
    options->weighting = options->weighting_conversion()->convert(weighting_string);
    
    // Get integral options
    options->integration_ordinates = input_node.get_child_value<int>("integration_ordinates");
    options->external_integral_calculation = input_node.get_attribute<bool>("external_integral_calculation",
                                                                            options->external_integral_calculation);
    options->perform_integration = false;
    
    if (options->external_integral_calculation)
    {
        options->adaptive_quadrature = input_node.get_attribute<bool>("adaptive_quadrature", options->adaptive_quadrature);
        if (options->adaptive_quadrature)
        {
            options->minimum_radius_ordinates = input_node.get_attribute<int>("minimum_radius_ordinates",
                                                                                 options->minimum_radius_ordinates);
            options->maximum_integration_ordinates = input_node.get_attribute<int>("maximum_integration_ordinates",
                                                                                 options->maximum_integration_ordinates);
        }
        
        options->dimensional_cells = input_node.get_child_vector<int>("dimensional_cells", dimension);
    }
    options->solid = solid_geometry_;
    meshless_factory.get_boundary_limits(dimension,
                                         boundary_surfaces_,
                                         options->limits);
    
    // Get Galerkin option
    options->identical_basis_functions = options->identical_basis_functions_conversion()->convert("false");
    
    // Get SUPG options
    options->include_supg = false;
    
    return options;
}

