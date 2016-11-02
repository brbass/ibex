#include "Boundary_Source_Parser.hh"

using namespace std;

Boundary_Source_Parser::
Boundary_Source_Parser(pugi::xml_node &input_file,
                       shared_ptr<Angular_Discretization> angular,
                       shared_ptr<Energy_Discretization> energy):
    Vector_Parser(input_file),
    angular_(angular),
    energy_(energy)
{
    pugi::xml_node sources_node = input_file.child("boundary_sources");
    
    int number_of_ordinates = angular_->number_of_ordinates();
    int number_of_groups = energy_->number_of_groups();

    int number_of_boundary_sources = XML_Functions::child_value<int>(sources_node, "number_of_boundary_sources");
    
    sources_.resize(number_of_boundary_sources);

    for (pugi::xml_node source_node = sources_node.child("boundary_source"); source_node; source_node = source_node.next_sibling("boundary_source"))
    {
        int a = XML_Functions::child_value<int>(source_node, "index");
        
        vector<double> alpha = XML_Functions::child_vector<double>(source_node, "alpha", number_of_groups);
        vector<double> isotropic_boundary_source = XML_Functions::child_vector<double>(source_node, "isotropic_source", number_of_groups);
        vector<double> boundary_source(number_of_groups * number_of_ordinates);
        
        for (int g = 0; g < number_of_groups; ++g)
        {
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                int k = g + number_of_groups * o;
                
                boundary_source[k] = isotropic_boundary_source[g];
            }
        }
        
        shared_ptr<Boundary_Source> source = make_shared<Boundary_Source>(a,
                                                                          angular_,
                                                                          energy_,
                                                                          boundary_source,
                                                                          alpha);
        
        sources_[a] = source;
    }
}
