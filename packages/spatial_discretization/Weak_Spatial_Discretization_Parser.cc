#include "Weak_Spatial_Discretization_Parser.hh"

Weak_Spatial_Discretization_Parser::
Weak_Spatial_Discretization_Parser()
{
}

shared_ptr<Weak_Spatial_Discretization> Weak_Spatial_Discretization_Parser::
parse_from_xml(XML_Node input_node) const
{
    // Get solid geometry


    // Get meshless function
    
    int number_of_points = input_node.get_child_value("number_of_points");
    
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
