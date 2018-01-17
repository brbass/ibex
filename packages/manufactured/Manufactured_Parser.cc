#include "Manufactured_Parser.hh"

#include <iostream>

#include "Angular_Discretization.hh"
#include "Cartesian_Plane.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "Manufactured_Constant_Cross_Sections.hh"
#include "Manufactured_Constant_Solution.hh"
#include "Manufactured_Solid_Geometry.hh"
#include "XML_Node.hh"

using namespace std;

Manufactured_Parser::
Manufactured_Parser(shared_ptr<Angular_Discretization> angular,
                    shared_ptr<Energy_Discretization> energy):
    angular_(angular),
    energy_(energy)
{
}

void Manufactured_Parser::
parse_from_xml(XML_Node input_node,
               shared_ptr<Manufactured_Solution> &solution,
               shared_ptr<Solid_Geometry> &solid,
               vector<shared_ptr<Cartesian_Plane> > &boundary_surfaces)
{
    solution = get_solution(input_node.get_child("solution"));
    shared_ptr<Manufactured_Cross_Sections> cross_sections = get_cross_sections(input_node.get_child("cross_sections"));
    solid = make_shared<Manufactured_Solid_Geometry>(angular_,
                                                     energy_,
                                                     solution,
                                                     cross_sections);
    boundary_surfaces = get_boundary_surfaces(input_node.get_child("boundary_surfaces"));
}

shared_ptr<Manufactured_Solution> Manufactured_Parser::
get_solution(XML_Node input_node)
{
    shared_ptr<Manufactured_Solution> solution;
    
    string solution_type = input_node.get_attribute<string>("type");
    if (solution_type == "constant")
    {
        int expected_size = (energy_->number_of_groups()
                             * angular_->number_of_moments());
        vector<double> data = input_node.get_child_vector<double>("values",
                                                                  expected_size);
        solution = make_shared<Manufactured_Constant_Solution>(angular_,
                                                               energy_,
                                                               data);
    }
    else
    {
        AssertMsg(false, "solution type (" + solution_type + ") not found");
    }
    
    return solution;
}

shared_ptr<Manufactured_Cross_Sections> Manufactured_Parser::
get_cross_sections(XML_Node input_node)
{
    shared_ptr<Manufactured_Cross_Sections> cross_sections;
    
    string type = input_node.get_attribute<string>("type");
    if (type == "constant")
    {
        int number_of_groups = energy_->number_of_groups();
        int number_of_scattering_moments = angular_->number_of_scattering_moments();
        vector<double> sigma_t = input_node.get_child_vector<double>("sigma_t",
                                                                     number_of_groups);
        vector<double> sigma_s = input_node.get_child_vector<double>("sigma_s",
                                                                     number_of_groups * number_of_groups * number_of_scattering_moments);
        cross_sections
            = make_shared<Manufactured_Constant_Cross_Sections>(angular_,
                                                                energy_,
                                                                sigma_t,
                                                                sigma_s);
    }
    else
    {
        AssertMsg(false, "cross section type (" + type + ") not found");
    }

    return cross_sections;
}

vector<shared_ptr<Cartesian_Plane> > Manufactured_Parser::
get_boundary_surfaces(XML_Node input_node)
{
    int dimension = angular_->dimension();
    vector<double> positions = input_node.get_child_vector<double>("positions",
                                                                   2 * dimension);
    vector<shared_ptr<Cartesian_Plane> > surfaces(2 * dimension);
    for (int d = 0; d < dimension; ++d)
    {
        for (int i = 0; i < 2; ++i)
        {
            double normal = i == 0 ? -1 : 1;
            int index = i + 2 * d;
            surfaces[index] = make_shared<Cartesian_Plane>(index,
                                                           dimension,
                                                           Surface::Surface_Type::BOUNDARY,
                                                           d,
                                                           positions[index],
                                                           normal);
        }
    }

    return surfaces;
}
