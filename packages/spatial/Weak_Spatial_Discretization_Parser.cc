#include "Weak_Spatial_Discretization_Parser.hh"

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
#include "RBF_Function.hh"
#include "RBF_Parser.hh"
#include "Solid_Geometry.hh"
#include "Truncated_Gaussian_RBF.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Factory.hh"
#include "Weight_Function.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;

Weak_Spatial_Discretization_Parser::
Weak_Spatial_Discretization_Parser(shared_ptr<Solid_Geometry> solid_geometry,
                                   vector<shared_ptr<Cartesian_Plane> > const &boundary_surfaces):
    solid_geometry_(solid_geometry),
    boundary_surfaces_(boundary_surfaces)
{
}

shared_ptr<Weak_Spatial_Discretization> Weak_Spatial_Discretization_Parser::
get_weak_discretization(XML_Node input_node) const
{
    // Initialize data
    string input_format = input_node.get_attribute<string>("input_format");
    
    // Get basis and weight functions that include all pertinant connectivity
    if (input_format == "full")
    {
        return get_full_discretization(input_node);
    }
    // Get basis and weight functions only given points
    else if (input_format == "galerkin_points")
    {
        return get_galerkin_points_discretization(input_node);
    }
    else if (input_format == "cartesian")
    {
        return get_cartesian_discretization(input_node);
    }
    else if (input_format == "legendre")
    {
        return get_legendre_discretization(input_node);
    }
    else
    {
        AssertMsg(false, "input format (" + input_format + ") not found");
        return shared_ptr<Weak_Spatial_Discretization>();
    }
}

shared_ptr<Weak_Spatial_Discretization> Weak_Spatial_Discretization_Parser::
get_full_discretization(XML_Node input_node) const
{
    Meshless_Function_Factory meshless_factory;
    
    // Initialize data
    int dimension = solid_geometry_->dimension();
    int number_of_points = input_node.get_child_value<int>("number_of_points");

    // Get discretization options
    shared_ptr<Weak_Spatial_Discretization_Options> weak_options
        = get_weak_options(input_node.get_child("options"));
    
    // Get dimensional moments
    shared_ptr<Dimensional_Moments> dimensional_moments
        = make_shared<Dimensional_Moments>(weak_options->include_supg,
                                           dimension);

    // Get basis and weight functions that include all pertinant connectivity
    vector<shared_ptr<Basis_Function> > basis_functions
        = get_basis_functions(input_node.get_child("basis_functions"),
                              number_of_points,
                              dimension);
    vector<shared_ptr<Weight_Function> > weight_functions
        = get_weight_functions(input_node.get_child("weight_functions"),
                               number_of_points,
                               dimension,
                               weak_options,
                               dimensional_moments,
                               basis_functions);
    
    // Create spatial discretization
    return make_shared<Weak_Spatial_Discretization>(basis_functions,
                                                    weight_functions,
                                                    dimensional_moments,
                                                    weak_options);
}

shared_ptr<Weak_Spatial_Discretization> Weak_Spatial_Discretization_Parser::
get_galerkin_points_discretization(XML_Node input_node) const
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
    Assert(weak_options->identical_basis_functions
           != Weak_Spatial_Discretization_Options::Identical_Basis_Functions::FALSE);
    weak_options->identical_basis_functions
        = Weak_Spatial_Discretization_Options::Identical_Basis_Functions::TRUE;
    
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
    
    // Get RBFs
    shared_ptr<RBF> rbf = rbf_parser.parse_from_xml(weights_node.get_child("meshless_function"));
    shared_ptr<Distance> distance = make_shared<Cartesian_Distance>(dimension);
    bool global_rbf = rbf->range() == RBF::Range::GLOBAL;
    
    // Get independent meshless functions
    vector<shared_ptr<Meshless_Function> > meshless_functions;
    meshless_factory.get_rbf_functions(number_of_points,
                                       radii,
                                       points,
                                       rbf,
                                       distance,
                                       meshless_functions);
    
    // Find neighbors
    vector<vector<int> > neighbors;
    vector<vector<double> > squared_distances;
    meshless_factory.get_neighbors(kd_tree,
                                   global_rbf,
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
    
    // Get MLS
    string meshless_type = weights_node.get_child("meshless_function").get_attribute<string>("type");
    if (meshless_type == "rbf")
    {
        // Do nothing
    }
    else if (meshless_type == "linear_mls")
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

    // Get basis functions
    vector<shared_ptr<Basis_Function> > basis_functions;
    weak_factory.get_basis_functions(number_of_points,
                                     meshless_functions,
                                     basis_functions);

    // Get weight functions
    vector<shared_ptr<Weight_Function> > weight_functions;
    weak_factory.get_weight_functions(number_of_points,
                                      weight_options,
                                      weak_options,
                                      dimensional_moments,
                                      neighbors,
                                      meshless_functions,
                                      basis_functions,
                                      weight_functions);
    
    // Create spatial discretization
    return make_shared<Weak_Spatial_Discretization>(basis_functions,
                                                    weight_functions,
                                                    dimensional_moments,
                                                    weak_options,
                                                    kd_tree);
}

shared_ptr<Weak_Spatial_Discretization> Weak_Spatial_Discretization_Parser::
get_cartesian_discretization(XML_Node input_node) const
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
    Assert(weak_options->identical_basis_functions
           != Weak_Spatial_Discretization_Options::Identical_Basis_Functions::FALSE);
    weak_options->identical_basis_functions
        = Weak_Spatial_Discretization_Options::Identical_Basis_Functions::TRUE;
    
    // Initialize data
    int dimension = solid_geometry_->dimension();
    XML_Node weights_node = input_node.get_child("weight_functions");
    
    // Get dimensional moments
    shared_ptr<Dimensional_Moments> dimensional_moments
        = make_shared<Dimensional_Moments>(weak_options->include_supg,
                                           dimension);
    
    // Get points
    vector<int> dimensional_points =
        input_node.get_child_vector<int>("dimensional_points",
                                         dimension);
    int number_of_points;
    vector<vector<double> > points;
    meshless_factory.get_cartesian_points(dimension,
                                          dimensional_points,
                                          weak_options->limits,
                                          number_of_points,
                                          points);
    
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
    bool global_rbf = rbf->range() == RBF::Range::GLOBAL;
    
    // Get independent meshless functions
    vector<shared_ptr<Meshless_Function> > meshless_functions;
    meshless_factory.get_rbf_functions(number_of_points,
                                       radii,
                                       points,
                                       rbf,
                                       distance,
                                       meshless_functions);
    
    // Find neighbors
    vector<vector<int> > neighbors;
    vector<vector<double> > squared_distances;
    meshless_factory.get_neighbors(kd_tree,
                                   global_rbf,
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
    
    // Get MLS
    string meshless_type = weights_node.get_child("meshless_function").get_attribute<string>("type");
    if (meshless_type == "rbf")
    {
        // Do nothing
    }
    else if (meshless_type == "linear_mls")
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

    // Get basis functions
    vector<shared_ptr<Basis_Function> > basis_functions;
    weak_factory.get_basis_functions(number_of_points,
                                     meshless_functions,
                                     basis_functions);

    // Get weight functions
    vector<shared_ptr<Weight_Function> > weight_functions;
    weak_factory.get_weight_functions(number_of_points,
                                      weight_options,
                                      weak_options,
                                      dimensional_moments,
                                      neighbors,
                                      meshless_functions,
                                      basis_functions,
                                      weight_functions);
    
    // Create spatial discretization
    return make_shared<Weak_Spatial_Discretization>(basis_functions,
                                                    weight_functions,
                                                    dimensional_moments,
                                                    weak_options,
                                                    kd_tree);
}

shared_ptr<Weak_Spatial_Discretization> Weak_Spatial_Discretization_Parser::
get_legendre_discretization(XML_Node input_node) const
{
    Meshless_Function_Factory meshless_factory;
    Weak_Spatial_Discretization_Factory weak_factory(solid_geometry_,
                                                     boundary_surfaces_);
    
    // Get options
    shared_ptr<Weight_Function_Options> weight_options
        = get_weight_options(input_node.get_child("options"));
    shared_ptr<Weak_Spatial_Discretization_Options> weak_options
        = get_weak_options(input_node.get_child("options"));
    Assert(weak_options->identical_basis_functions
           != Weak_Spatial_Discretization_Options::Identical_Basis_Functions::FALSE);
    weak_options->identical_basis_functions
        = Weak_Spatial_Discretization_Options::Identical_Basis_Functions::TRUE;
    
    // Initialize data
    int dimension = solid_geometry_->dimension();
    
    // Get dimensional moments
    shared_ptr<Dimensional_Moments> dimensional_moments
        = make_shared<Dimensional_Moments>(weak_options->include_supg,
                                           dimension);
    
    // Get Legendre functions
    vector<int> order
        = input_node.get_child_vector<int>("order",
                                           dimension);
    int number_of_points;
    vector<shared_ptr<Meshless_Function> > meshless_functions;
    meshless_factory.get_legendre_functions(dimension,
                                            order,
                                            weak_options->limits,
                                            number_of_points,
                                            meshless_functions);

    // Set all functions to be neighbors of one another
    vector<int> local_neighbors(number_of_points);
    for (int i = 0; i < number_of_points; ++i)
    {
        local_neighbors[i] = i;
    }
    vector<vector<int> > neighbors(number_of_points, local_neighbors);

    // Set each function to come first in its own list
    for (int i = 0; i < number_of_points; ++i)
    {
        neighbors[i][0] = i;
        neighbors[i][i] = 0;
    }
    
    // Get MLS
    string meshless_type = input_node.get_child("meshless_function").get_attribute<string>("type");
    if (meshless_type == "legendre")
    {
        // Do nothing
    }
    else if (meshless_type == "linear_mls")
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
    
    // Get basis functions
    vector<shared_ptr<Basis_Function> > basis_functions;
    weak_factory.get_basis_functions(number_of_points,
                                     meshless_functions,
                                     basis_functions);

    // Get weight functions
    vector<shared_ptr<Weight_Function> > weight_functions;
    weak_factory.get_weight_functions(number_of_points,
                                      weight_options,
                                      weak_options,
                                      dimensional_moments,
                                      neighbors,
                                      meshless_functions,
                                      basis_functions,
                                      weight_functions);
    
    // Create spatial discretization
    return make_shared<Weak_Spatial_Discretization>(basis_functions,
                                                    weight_functions,
                                                    dimensional_moments,
                                                    weak_options);
}

vector<shared_ptr<Meshless_Function> > Weak_Spatial_Discretization_Parser::
get_rbf_functions(XML_Node input_node,
                  int number_of_points,
                  int dimension,
                  string prefix) const
{
    // Get RBF and distance information
    RBF_Parser rbf_parser;
    shared_ptr<RBF> rbf = rbf_parser.parse_from_xml(input_node.get_child("meshless_function"));
    shared_ptr<Distance> distance = make_shared<Cartesian_Distance>(dimension);
    
    // Get meshless functions
    vector<shared_ptr<Meshless_Function> > functions(number_of_points);
    for (XML_Node node = input_node.get_child(prefix);
         node;
         node = node.get_sibling(prefix,
                                 false))
    {
        int index = node.get_attribute<int>("index");
        Assert(index < number_of_points);
        double radius = node.get_child_value<double>("radius");
        vector<double> position = node.get_child_vector<double>("position",
                                                                dimension);
        double shape = rbf->radius() / radius;
        functions[index] = make_shared<RBF_Function>(index,
                                                     shape,
                                                     position,
                                                     rbf,
                                                     distance);
    }

    return functions;
}

vector<shared_ptr<Meshless_Function> > Weak_Spatial_Discretization_Parser::
get_mls_functions(XML_Node input_node,
                  int order,
                  int number_of_points,
                  int dimension,
                  string prefix) const
{
    // Get weighting functions for MLS
    vector<shared_ptr<Meshless_Function> > rbf_functions
        = get_rbf_functions(input_node,
                            number_of_points,
                            dimension,
                            prefix);

    
    // Get MLS functions
    vector<shared_ptr<Meshless_Function> > functions(number_of_points);
    for (XML_Node node = input_node.get_child(prefix);
         node;
         node = node.get_sibling(prefix,
                                 false))
    {
        int index = node.get_attribute<int>("index");
        
        int number_of_neighbors
            = node.get_child_value<int>("number_of_" + prefix + "_neighbors");
        vector<int> neighbor_indices
            = node.get_child_vector<int>(prefix + "_neighbors",
                                                number_of_neighbors);
        
        vector<shared_ptr<Meshless_Function> > local_functions(number_of_neighbors);
        for (int i = 0; i < number_of_neighbors; ++i)
        {
            local_functions[i] = rbf_functions[neighbor_indices[i]];
        }

        switch (order)
        {
        case 1:
            functions[index] = make_shared<Linear_MLS_Function>(local_functions);
            break;
        case 2:
            functions[index] = make_shared<Quadratic_MLS_Function>(local_functions);
            break;
        default:
            AssertMsg(false, "order for MLS function not found");
        }
    }

    return functions;
}

vector<shared_ptr<Meshless_Function> > Weak_Spatial_Discretization_Parser::
get_meshless_functions(XML_Node input_node,
                       int number_of_points,
                       int dimension,
                       string prefix) const
{
    XML_Node meshless_node = input_node.get_child("meshless_function");
    string meshless_type = meshless_node.get_attribute<string>("type");
    if (meshless_type == "rbf")
    {
        return get_rbf_functions(input_node,
                                 number_of_points,
                                 dimension,
                                 prefix);
    }
    else if (meshless_type == "linear_mls")
    {
        return get_mls_functions(input_node,
                                 1,
                                 number_of_points,
                                 dimension,
                                 prefix);
    }
    else if (meshless_type == "quadratic_mls")
    {
        return get_mls_functions(input_node,
                                 2,
                                 number_of_points,
                                 dimension,
                                 prefix);
    }
    else
    {
        AssertMsg(false, "meshless type (" + meshless_type + ") not found");
        return vector<shared_ptr<Meshless_Function> >();
    }
}

vector<shared_ptr<Basis_Function> >  Weak_Spatial_Discretization_Parser::
get_basis_functions(XML_Node input_node,
                    int number_of_points,
                    int dimension) const
{
    // Get meshless functions
    vector<shared_ptr<Meshless_Function> > meshless_functions
        = get_meshless_functions(input_node,
                                 number_of_points,
                                 dimension,
                                 "basis");

    // Get basis functions
    vector<shared_ptr<Basis_Function> > basis_functions(number_of_points);
    for (XML_Node node = input_node.get_child("basis");
         node;
         node = node.get_sibling("basis",
                                 false))
    {
        int index = node.get_attribute<int>("index");

        // Get local boundary surfaces
        vector<shared_ptr<Cartesian_Plane> > boundary_surfaces
            = get_boundary_surfaces(meshless_functions[index]);

        // Create basis function
        basis_functions[index] = make_shared<Basis_Function>(index,
                                                             dimension,
                                                             meshless_functions[index],
                                                             boundary_surfaces);
    }

    return basis_functions;
}

shared_ptr<Weight_Function_Options> Weak_Spatial_Discretization_Parser::
get_weight_options(XML_Node input_node) const
{
    shared_ptr<Weight_Function_Options> options
        = make_shared<Weight_Function_Options>();
    
    options->tau_const
        = input_node.get_child_value<double>("tau",
                                             options->tau_const);
    options->output_material
        = input_node.get_attribute<bool>("output_material",
                                         options->output_material);
    options->output_integrals
        = input_node.get_attribute<bool>("output_integrals",
                                         options->output_integrals);

    return options;
}

shared_ptr<Weak_Spatial_Discretization_Options> Weak_Spatial_Discretization_Parser::
get_weak_options(XML_Node input_node) const
{
    int dimension = solid_geometry_->dimension();
    shared_ptr<Weak_Spatial_Discretization_Options> options
        = make_shared<Weak_Spatial_Discretization_Options>();
    options->discretization = Weak_Spatial_Discretization_Options::Discretization::WEAK;
    Meshless_Function_Factory meshless_factory;
    
    // Get weighting method
    string weighting_string = input_node.get_attribute<string>("weighting");
    options->weighting = options->weighting_conversion()->convert(weighting_string);
    if (options->weighting == Weak_Spatial_Discretization_Options::Weighting::FLUX)
    {
        string coeff_filename = input_node.get_attribute<string>("flux_file");
        vector<string> coeff_path = input_node.get_attribute_vector<string>("flux_path");
        XML_Document coeff_doc(coeff_filename);
        XML_Node coeff_node = coeff_doc.get_child(coeff_path[0]);
        for (int i = 1; i < coeff_path.size(); ++i)
        {
            coeff_node = coeff_node.get_child(coeff_path[i]);
        }
        options->flux_coefficients = coeff_node.get_vector<double>();
        options->scalar_flux_fraction = input_node.get_attribute<double>("scalar_flux_fraction");
    }
    
    // Get integral options
    options->integration_ordinates = input_node.get_child_value<int>("integration_ordinates");
    options->external_integral_calculation = input_node.get_attribute<bool>("external_integral_calculation",
                                                                            options->external_integral_calculation);
    options->perform_integration = input_node.get_attribute<bool>("perform_integration",
options->perform_integration);
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
    string identical_string = input_node.get_attribute<string>("identical_basis_functions");
    options->identical_basis_functions = options->identical_basis_functions_conversion()->convert(identical_string);
    
    // Get SUPG options
    options->include_supg = input_node.get_attribute<bool>("supg");
    
    // Get tau scaling
    string tau_string = input_node.get_attribute<string>("tau_scaling");
    options->tau_scaling = options->tau_scaling_conversion()->convert(tau_string);
    
    return options;
}

vector<shared_ptr<Weight_Function> > Weak_Spatial_Discretization_Parser::
get_weight_functions(XML_Node input_node,
                     int number_of_points,
                     int dimension,
                     shared_ptr<Weak_Spatial_Discretization_Options> weak_options,
                     shared_ptr<Dimensional_Moments> dimensional_moments,
                     vector<shared_ptr<Basis_Function> > const &basis_functions) const
{
    // Get global weight function information
    shared_ptr<Weight_Function_Options> weight_options = get_weight_options(input_node);
    
    // Get meshless functions
    vector<shared_ptr<Meshless_Function> > meshless_functions
        = get_meshless_functions(input_node,
                                 number_of_points,
                                 dimension,
                                 "weight");
    
    // Get weight functions
    vector<shared_ptr<Weight_Function> > weight_functions(number_of_points);
    for (XML_Node node = input_node.get_child("weight");
         node;
         node = node.get_sibling("weight",
                                 false))
    {
        int index = node.get_attribute<int>("index");
        
        // Get local basis functions
        int number_of_basis_neighbors
            = node.get_child_value<int>("number_of_basis_neighbors");
        vector<int> basis_neighbor_indices
            = node.get_child_vector<int>("basis_neighbors",
                                         number_of_basis_neighbors);
        vector<shared_ptr<Basis_Function> > local_basis_functions(number_of_basis_neighbors);
        for (int i = 0; i < number_of_basis_neighbors; ++i)
        {
            local_basis_functions[i] = basis_functions[basis_neighbor_indices[i]];
        }

        // Get local boundary surfaces
        vector<shared_ptr<Cartesian_Plane> > boundary_surfaces
            = get_boundary_surfaces(meshless_functions[index]);

        // Create weight function
        weight_functions[index]
            = make_shared<Weight_Function>(index,
                                           dimension,
                                           weight_options,
                                           weak_options,
                                           meshless_functions[index],
                                           local_basis_functions,
                                           dimensional_moments,
                                           solid_geometry_,
                                           boundary_surfaces);
    }

    return weight_functions;
}

vector<shared_ptr<Cartesian_Plane> > Weak_Spatial_Discretization_Parser::
get_boundary_surfaces(shared_ptr<Meshless_Function> function) const
{
    vector<shared_ptr<Cartesian_Plane> > surfaces;
    
    for (shared_ptr<Cartesian_Plane> surface : boundary_surfaces_)
    {
        double surface_position = surface->position();
        int surface_dim = surface->surface_dimension();
        vector<double> function_position = function->position();
        double function_radius = function->radius();

        if (std::abs(function_position[surface_dim] - surface_position) <= function_radius)
        {
            surfaces.push_back(surface);
        }
    }
    
    return surfaces;
}

