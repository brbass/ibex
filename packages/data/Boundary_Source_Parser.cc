#include "Boundary_Source_Parser.hh"

#include "Angular_Discretization.hh"
#include "Boundary_Source.hh"
#include "Energy_Discretization.hh"
#include "XML_Node.hh"

using namespace std;

Boundary_Source_Parser::
Boundary_Source_Parser(shared_ptr<Angular_Discretization> angular,
                       shared_ptr<Energy_Discretization> energy):
    angular_(angular),
    energy_(energy)
{
}

vector<shared_ptr<Boundary_Source> > Boundary_Source_Parser::
parse_from_xml(XML_Node input_node)
{
    int number_of_ordinates = angular_->number_of_ordinates();
    int number_of_groups = energy_->number_of_groups();

    int number_of_boundary_sources = input_node.get_child_value<int>("number_of_boundary_sources");
    
    vector<shared_ptr<Boundary_Source> > sources(number_of_boundary_sources);

    int checksum = 0;
    for (XML_Node source_node = input_node.get_child("boundary_source"); source_node; source_node = source_node.get_sibling("boundary_source"))
    {
        int a = source_node.get_child_value<int>("index");
        
        vector<double> alpha = source_node.get_child_vector<double>("alpha", number_of_groups);
        vector<double> isotropic_boundary_source = source_node.get_child_vector<double>("isotropic_source", number_of_groups);
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
        
        sources[a] = source;
        
        checksum += 1;
    } // boundary sources

    int checksum_expected = number_of_boundary_sources * (number_of_boundary_sources - 1) / 2;
    AssertMsg(checksum == checksum_expected, "Boundary source indexing incorrect");

    return sources;
}
