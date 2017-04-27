#include "Weak_Spatial_Discretization_Parser.hh"

#include "Basis_Function.hh"
#include "Cartesian_Distance.hh"
#include "Cartesian_Plane.hh"
#include "Compact_Gaussian_RBF.hh"
#include "Dimensional_Moments.hh"
#include "Linear_MLS_Function.hh"
#include "Meshless_Function.hh"
#include "RBF_Function.hh"
#include "RBF_Parser.hh"
#include "Solid_Geometry.hh"
#include "Truncated_Gaussian_RBF.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weight_Function.hh"
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
    // Get data
    int dimension = solid_geometry_->dimension();
    int number_of_points = input_node.get_child_value<int>("number_of_points");
    vector<shared_ptr<Basis_Function> > basis_functions
        = get_basis_functions(input_node.get_child("basis_functions"),
                              number_of_points,
                              dimension);
    vector<shared_ptr<Weight_Function> > weight_functions
        = get_weight_functions(input_node.get_child("weight_functions"),
                               number_of_points,
                               dimension,
                               basis_functions);
    bool supg = weight_functions[0]->options().include_supg;
    shared_ptr<Dimensional_Moments> dimensional_moments
        = make_shared<Dimensional_Moments>(supg,
                                           dimension);
    Weak_Spatial_Discretization::Options options;
    if (weight_functions[0]->options().external_integral_calculation)
    {
        AssertMsg(false, "not implemented for parser");
    }
    return make_shared<Weak_Spatial_Discretization>(basis_functions,
                                                    weight_functions,
                                                    dimensional_moments,
                                                    options);
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
        double radius = node.get_child_value<double>("radius");
        vector<double> position = node.get_child_vector<double>("position",
                                                                dimension);
        double shape = rbf->radius() / radius;
        functions[index] = make_shared<RBF_Function>(shape,
                                                     position,
                                                     rbf,
                                                     distance);
    }

    return functions;
}

vector<shared_ptr<Meshless_Function> > Weak_Spatial_Discretization_Parser::
get_mls_functions(XML_Node input_node,
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
        
        functions[index] = make_shared<Linear_MLS_Function>(local_functions);
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

vector<shared_ptr<Weight_Function> > Weak_Spatial_Discretization_Parser::
get_weight_functions(XML_Node input_node,
                     int number_of_points,
                     int dimension,
                     vector<shared_ptr<Basis_Function> > const &basis_functions) const
{
    // Get global weight function information
    Weight_Function::Options material_options;
    {
        int integration_ordinates = input_node.get_child_value<int>("integration_ordinates");
        material_options.integration_ordinates = integration_ordinates;
        XML_Node material_node = input_node.get_child("material_options");
        string weighting = material_node.get_attribute<string>("weighting");
        if (weighting == "point")
        {
            material_options.weighting = Weight_Function::Options::Weighting::POINT;
        }
        else if (weighting == "weight")
        {
            material_options.weighting = Weight_Function::Options::Weighting::WEIGHT;
        }
        else if (weighting == "flux")
        {
            material_options.weighting = Weight_Function::Options::Weighting::FLUX;
        }
        else
        {
            AssertMsg(false, "weighting option (" + weighting + ") not found");
        }
        
        string output = material_node.get_attribute<string>("output");
        if (output == "standard")
        {
            material_options.output = Weight_Function::Options::Output::STANDARD;
        }
        else if (output == "supg")
        {
            material_options.output = Weight_Function::Options::Output::SUPG;
        }
        else
        {
            AssertMsg(false, "output option (" + output + ") not found");
        }

        material_options.tau_const = material_node.get_attribute<double>("tau", 1.0);
    }
    
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
                                           material_options,
                                           meshless_functions[index],
                                           local_basis_functions,
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

        if (abs(function_position[surface_dim] - surface_position) <= function_radius)
        {
            surfaces.push_back(surface);
        }
    }

    return surfaces;
}

