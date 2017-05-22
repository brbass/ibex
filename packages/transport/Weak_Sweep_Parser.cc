#include "Weak_Sweep_Parser.hh"

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Transport_Discretization.hh"
#include "Weak_RBF_Sweep.hh"
#include "Weak_Spatial_Discretization.hh"
#include "XML_Node.hh"

using namespace std;

Weak_Sweep_Parser::
Weak_Sweep_Parser(shared_ptr<Weak_Spatial_Discretization> spatial,
                  shared_ptr<Angular_Discretization> angular,
                  shared_ptr<Energy_Discretization> energy,
                  shared_ptr<Transport_Discretization> transport):
    spatial_(spatial),
    angular_(angular),
    energy_(energy),
    transport_(transport)
{
}

shared_ptr<Weak_RBF_Sweep> Weak_Sweep_Parser::
get_weak_rbf_sweep(XML_Node input_node) const
{
    Weak_RBF_Sweep::Options options;
    
    string solver = input_node.get_attribute<string>("solver", "amesos");

    if (solver == "amesos")
    {
        options.solver = Weak_RBF_Sweep::Options::Solver::AMESOS;
    }
    else if (solver == "aztec")
    {
        options.solver = Weak_RBF_Sweep::Options::Solver::AZTEC;
    }
    else if (solver == "aztec_ilut")
    {
        options.solver = Weak_RBF_Sweep::Options::Solver::AZTEC_ILUT;
    }
    else
    {
        AssertMsg(false, "solver (" + solver + ") not found");
    }
    
    return make_shared<Weak_RBF_Sweep>(options,
                                       spatial_,
                                       angular_,
                                       energy_,
                                       transport_);
}
