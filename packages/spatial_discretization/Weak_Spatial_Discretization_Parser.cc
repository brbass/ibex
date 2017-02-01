#include "Weak_Spatial_Discretization_Parser.hh"

#include "Boundary_Source.hh"
#include "Compact_Gaussian_RBF.hh"
#include "Material.hh"
#include "Meshless_Function.hh"
#include "Solid_Geometry.hh"
#include "Truncated_Gaussian_RBF.hh"

Weak_Spatial_Discretization_Parser::
Weak_Spatial_Discretization_Parser(vector<shared_ptr<Material> > const &materials,
                                   vector<shared_ptr<Boundary_Source> > const &boundary_sources,
                                   shared_ptr<Solid_Geometry> solid_geometry)
{
}

shared_ptr<Weak_Spatial_Discretization> Weak_Spatial_Discretization_Parser::
parse_from_xml(XML_Node input_node) const
{
    int number_of_points = input_node.get_child_value("number_of_points");

    // Get meshless function information
    XML_Node meshless_node = input_node.get_child("meshless");
    string meshless_type = meshless_node.get_attribute<string>("type");
    
    string rbf_type = meshless_node.get_attribute<string>("function");
    shared_ptr<RBF> rbf;
    
    if (rbf_type == "compact_gaussian")
    {
        double radius = meshless_node.get_child_value<double>("radius",
                                                              5.0);
        rbf = make_shared<Compact_Gaussian_RBF>(radius);
    }
    else if (rbf_type == "truncated_gaussian")
    {
        double radius = meshless_node.get_child_value<double>("radius",
                                                              5.0);
        rbf = make_shared<Truncated_Gaussian_RBF>(radius);
    }
    else if (rbf_type == "wendland30")
    {
        rbf = make_shared<Wendland_RBF>(0);
    }
    else if (rbf_type == "wendland31")
    {
        rbf = make_shared<Wendland_RBF>(1);
    }
    else if (rbf_type == "wendland32")
    {
        rbf = make_shared<Wendland_RBF>(2);
    }
    else if (rbf_type == "wendland33")
    {
        rbf = make_shared<Wendland_RBF>(3);
    }
    else
    {
        AssertMsg(false, "basis_type \"" + rbf_type + "\" not found");
    }

    // Get meshless functions
    XML_Node bases_node = input_node.get_child("basis_functions");
    vector<shared_ptr<Meshless_Function> > meshless_functions(number_of_points);
    
    
    
    // Parse basis functions
    XML_Node bases_node = input_node.get_child("basis_functions");
    vector<shared_ptr<Basis_Function> > basis_functions(number_of_points);
    for (XML_Node basis_node = input_node.get_child("basis");
         basis_node;
         basis_node = basis_node.get_sibling("basis",
                                             false))
    {
        int index = basis_node.get_attribute<int>("index");
        double max_distance = basis_node.get_child_value("max_distance");
        vector<double> position = basis_node.get_child_vector("position");
        shared_ptr<Meshless_Function> meshless_function
            = get_meshless_function();
        
        basis_functions[index]
            = make_shared<Basis_Function>(index,
                                          dimension,
                                          meshless_function,
                                          boundary_surfaces);
    }

    // Parse weight functions
    XML_Node weights_node = input_node.get_child("weight_functions");
    vector<shared_ptr<Weight_Function> > weight_functions(number_of_points);
    for (XML_Node weight_node = input_node.get_child("weight");
         weight_node;
         weight_node = weight_node.get_sibling("weight",
                                               false))
    {
        int index = weight_node.get_attribute<int>("index");
        double max_distance = weight_node.get_child_value("max_distance");
        vector<double> position = basis_node.get_child_vector("position");
        vector<int> neighbors = weight_node.get_child_vector("neighbors");
        
    }
    
    
    
}

shared_ptr<Meshless_Function> Weak_Spatial_Discretization_Parser::
get_meshless_function() const
{
}
