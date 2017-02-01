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
parse_from_xml(XML_Node input_node) const
{
    int number_of_points = input_node.get_child_value("number_of_points");
    int dimension = solid_geometry->dimension();
    
    XML_Node meshless_node = input_node.get_child("meshless");
    RBF_Parser rbf_parser;
    shared_ptr<RBF> rbf = rbf_parser->parse_from_xml(meshless_node);
    shared_ptr<Distance> distance = make_shared<Cartesian_Distance>(dimension);
    
    // Get meshless functions
    XML_Node bases_node = input_node.get_child("basis_functions");
    vector<shared_ptr<Meshless_Function> > meshless_basis(number_of_points);
    for (XML_Node basis_node = input_node.get_child("basis");
         basis_node;
         basis_node = basis_node.get_sibling("basis",
                                             false))
    {
        int index = basis_node.get_attribute<int>("index");
        double radius = basis_node.get_child_value("radius");
        vector<double> position = basis_node.get_child_vector("position");
        double shape = rbf->radius() / radius;
        meshless_basis[index] = make_shared<RBF_Function>(shape,
                                                          position,
                                                          rbf,
                                                          distance);
    }
    
    XML_Node weights_node = input_node.get_child("weight_functions");
    vector<shared_ptr<Meshless_Function> > meshless_weight(number_of_points);
    for (XML_Node weight_node = input_node.get_child("weight");
         weight_node;
         weight_node = weight_node.get_sibling("weight",
                                             false))
    {
        int index = weight_node.get_attribute<int>("index");
        double radius = weight_node.get_child_value("radius");
        vector<double> position = weight_node.get_child_vector("position");
        double shape = rbf->radius() / radius;
        meshless_weight[index] = make_shared<RBF_Function>(shape,
                                                          position,
                                                          rbf,
                                                          distance);
    }

    string meshless_type = meshless_node.get_attribute<string>("type");
    
    if (meshless_type == "rbf")
    {
        
    }
    else if (meshless_type == "linear_mls")
    {
        
    }
    else
    {
        
    }
}

shared_ptr<Meshless_Function> Weak_Spatial_Discretization_Parser::
get_meshless_function() const
{
}
