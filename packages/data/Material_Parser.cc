#include "Material_Parser.hh"

#include <iostream>

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"

using namespace std;

Material_Parser::
Material_Parser(pugi::xml_node &input_file,
                shared_ptr<Angular_Discretization> angular,
                shared_ptr<Energy_Discretization> energy):
    Vector_Parser(input_file),
    angular_(angular),
    energy_(energy)
{
    pugi::xml_node materials_node = input_file.child("materials");

    int number_of_materials = XML_Functions::child_value<int>(materials_node, "number_of_materials");
    
    int number_of_moments = angular_->number_of_scattering_moments();
    int number_of_groups = energy_->number_of_groups();
    
    // parse the data for each material

    materials_.resize(number_of_materials);
    
    for (pugi::xml_node material_node = materials_node.child("material"); material_node; material_node = material_node.next_sibling("material"))
    {
        int a = XML_Functions::child_value<int>(material_node, "index");
        
        vector<double> sigma_t = XML_Functions::child_vector<double>(material_node, "sigma_t", number_of_groups);
        vector<double> sigma_s = XML_Functions::child_vector<double>(material_node, "sigma_s", number_of_groups * number_of_groups * number_of_moments);
        vector<double> nu = XML_Functions::child_vector<double>(material_node, "nu", number_of_groups);
        vector<double> sigma_f = XML_Functions::child_vector<double>(material_node, "sigma_f", number_of_groups);
        vector<double> chi = XML_Functions::child_vector<double>(material_node, "chi", number_of_groups);
        vector<double> internal_source = XML_Functions::child_vector<double>(material_node, "internal_source", number_of_groups);
        
        shared_ptr<Material> material = make_shared<Material>(a,
                                                              angular_,
                                                              energy_,
                                                              sigma_t,
                                                              sigma_s,
                                                              nu,
                                                              sigma_f,
                                                              chi,
                                                              internal_source);

        materials_[a] = material;
    } // materials
} // constructor
