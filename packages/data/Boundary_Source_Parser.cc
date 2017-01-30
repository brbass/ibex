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
    for (XML_Node source_node = input_node.get_child("boundary_source");
         source_node;
         source_node = source_node.get_sibling("boundary_source",
                                               false))
    {
        int a = source_node.get_attribute<int>("index");
        
        vector<double> alpha = source_node.get_child_vector<double>("alpha", number_of_groups);
        vector<double> isotropic_boundary_source = source_node.get_child_vector<double>("isotropic_source", number_of_groups);
        
        Boundary_Source::Dependencies dependencies;
        dependencies.angular = Boundary_Source::Dependencies::Angular::ISOTROPIC;
        shared_ptr<Boundary_Source> source = make_shared<Boundary_Source>(a,
                                                                          dependencies,
                                                                          angular_,
                                                                          energy_,
                                                                          isotropic_boundary_source,
                                                                          alpha);
        
        sources[a] = source;
        
        checksum += a;
    } // boundary sources

    int checksum_expected = number_of_boundary_sources * (number_of_boundary_sources - 1) / 2;
    AssertMsg(checksum == checksum_expected, "Boundary source indexing incorrect");

    return sources;
}

shared_ptr<Boundary_Source> Boundary_Source_Parser::
get_vacuum_boundary()
{
    int number_of_ordinates = angular_->number_of_ordinates();
    int number_of_groups = energy_->number_of_groups();

    int index = -1; // add enum if this function is used

    vector<double> alpha(number_of_groups, 0);
    vector<double> boundary_source(number_of_groups * number_of_ordinates, 0);
    Boundary_Source::Dependencies dependencies;
    return make_shared<Boundary_Source>(index,
                                        dependencies,
                                        angular_,
                                        energy_,
                                        boundary_source,
                                        alpha);
}

shared_ptr<Boundary_Source> Boundary_Source_Parser::
get_reflective_boundary()
{
    int number_of_ordinates = angular_->number_of_ordinates();
    int number_of_groups = energy_->number_of_groups();

    int index = -2; // add enum if this function is used
    
    vector<double> alpha(number_of_groups, 1);
    vector<double> boundary_source(number_of_groups * number_of_ordinates, 0);
    
    Boundary_Source::Dependencies dependencies;
    return make_shared<Boundary_Source>(index,
                                        dependencies,
                                        angular_,
                                        energy_,
                                        boundary_source,
                                        alpha);
}
