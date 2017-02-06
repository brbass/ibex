#include "Weak_Spatial_Discretization_Parser.hh"

#include "Basis_Function.hh"
#include "Boundary_Source.hh"
#include "Compact_Gaussian_RBF.hh"
#include "Material.hh"
#include "Meshless_Function.hh"
#include "RBF_Parser.hh"
#include "Solid_Geometry.hh"
#include "Truncated_Gaussian_RBF.hh"
#include "Weight_Function.hh"

Weak_Spatial_Discretization_Parser::
Weak_Spatial_Discretization_Parser(shared_ptr<Solid_Geometry> solid_geometry,
                                   vector<shared_ptr<Cartesian_Plane> > boundary_surfaces):
    solid_geometry_(solid_geometry),
    boundary_surfaces_(boundary_surfaces)
{
}

shared_ptr<Weak_Spatial_Discretization> Weak_Spatial_Discretization_Parser::
get_weak_discretization(XML_Node input_node) const
{
    int dimension = solid_geometry->dimension();
    int number_of_points = input_node.get_child_value<int>("number_of_points");
    XML_Node basis_node = input_node.child("basis_functions");
    vector<shared_ptr<Basis_Function> > basis_functions
        = get_basis_functions(input_node,
                              number_of_points,
                              dimension);
    vector<shared_ptr<Weight_Function> > weight_functions
        = get_weight_functions(input_node,
                               number_of_points,
                               dimension,
                               basis_functions);

    return make_shared<Weak_Spatial_Discretization>(basis_functions,
                                                    weight_functions);
}

vector<shared_ptr<RBF_Function> > Weak_Spatial_Discretization_Parser::
get_rbf_functions(XML_Node input_node,
                  int number_of_points,
                  int dimension,
                  string prefix)
{
    // Get RBF and distance information
    RBF_Parser rbf_parser;
    shared_ptr<RBF> rbf = rbf_parser->parse_from_xml(input_node);
    shared_ptr<Distance> distance = make_shared<Cartesian_Distance>(dimension);
    
    // Get meshless functions
    vector<shared_ptr<RBF_Function> > functions(number_of_points);
    for (XML_Node node = input_node.get_child(prefix);
         node;
         node = node.get_sibling(prefix,
                                 false))
    {
        int index = node.get_attribute<int>("index");
        double radius = node.get_child_value("radius");
        vector<double> position = node.get_child_vector("position");
        double shape = rbf->radius() / radius;
        functions[index] = make_shared<RBF_Function>(shape,
                                                     position,
                                                     rbf,
                                                     distance);
    }

    return functions;
}

vector<shared_ptr<Linear_MLS_Function> > Weak_Spatial_Discretization_Parser::
get_mls_functions(XML_Node input_node,
                  int number_of_points,
                  int dimension,
                  string prefix)
{
    // Get weighting functions for MLS
    vector<shared_ptr<Meshless_Function> > rbf_functions
        = get_rbf_functions(input_node,
                            number_of_points,
                            dimension,
                            prefix);

    
    // Get MLS functions
    vector<shared_ptr<Linear_MLS_Function> > functions(number_of_points);
    for (XML_Node node = input_node.get_child(prefix);
         node;
         node = node.get_sibling(prefix,
                                 false))
    {
        int index = node.get_attribute<int>("index");
        
        int number_of_neighbors
            = weight_node.get_child_value("number_of_" + prefix + "_neighbors");
        vector<int> neighbor_indices
            = weight_node.get_child_vector<int>(prefix + "_neighbors",
                                                number_of_neighbors);
        
        vector<shared_ptr<Meshless_Function> > local_functions(number_of_basis_neighbors);
        for (int i = 0; i < number_of_basis_neighbors; ++i)
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
                       string prefix)
{
    string meshless_type = input_node.get_child_value<string>("meshless_type");
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
                    int dimension)
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
                     vector<shared_ptr<Basis_Function> > const &basis_functions)
{
    // Get global weight function information
    int integration_ordinates = input_node.get_child_value<int>("integration_ordinates");
    Weight_Function::Material_Options material_options;
    {
        XML_Node material_node = input_node.get_child("material_options");
        string weighting = material_node.get_child_value<string>("weighting",
                                                                 "weight");
        if (weighting == "point")
        {
            material_options.weighting = Weight_Function::Weighting::POINT;
        }
        else if (weighting == "weight")
        {
            material_options.weighting = Weight_Function::Weighting::WEIGHT;
        }
        else if (weighting == "weight")
        {
            material_options.weighting = Weight_Function::Weighting::FLUX;
        }
        else
        {
            AssertMsg(false, "weighting option (" + weighting + ") not found");
        }
        
        string output = material_node.get_child_value<string>("output",
                                                              "standard");
        if (output == "standard")
        {
            material_options.output = Weight_Function::Weighting::STANDARD;
        }
        else if (output == "standard")
        {
            material_options.output = Weight_Function::Weighting::SUPG;
        }
        else
        {
            AssertMsg(false, "output option (" + output + ") not found");
        }
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
            = weight_node.get_child_value("number_of_basis_neighbors");
        vector<int> basis_neighbor_indices
            = weight_node.get_child_vector<int>("basis_neighbors",
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
                                           integration_ordinates,
                                           material_options,
                                           meshless_weight[index],
                                           local_basis_functions,
                                           solid_geometry_,
                                           boundary_surfaces);
    }

    return weight_functions;
}

vector<shared_ptr<Cartesian_Plane> > Weak_Spatial_Discretization_Parser::
get_boundary_surfaces(shared_ptr<Meshless_Function> function)
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
