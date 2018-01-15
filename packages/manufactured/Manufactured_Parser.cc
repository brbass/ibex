#include "Manufactured_Parser.hh"

#include "Angular_Discretization.hh"
#include "Cartesian_Plane.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "Manufactured_Solution.hh"
#include "Solid_Geometry.hh"
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
    AssertMsg(false, "not yet implemented");
}
               
