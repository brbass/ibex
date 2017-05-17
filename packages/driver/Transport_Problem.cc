#include "Transport_Problem.hh"

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source_Parser.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Constructive_Solid_Geometry_Parser.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Krylov_Eigenvalue.hh"
#include "Krylov_Steady_State.hh"
#include "Material_Parser.hh"
#include "Solver.hh"
#include "Solver_Parser.hh"
#include "Source_Iteration.hh"
#include "Transport_Discretization.hh"
#include "Weak_RBF_Sweep.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "Weak_Sweep_Parser.hh"

using namespace std;

Transport_Problem::
Transport_Problem(XML_Node input_node,
                  XML_Node output_node):
    input_node_(input_node),
    output_node_(output_node)
{
}

void Transport_Problem::
solve()
{
    XML_Node problem_node = input_node_.get_child("problem");
    string type = problem_node.get_attribute<string>("type");
    
    if (type == "eigenvalue")
    {
        solve_eigenvalue();
    }
    else if (type == "steady_state")
    {
        solve_steady_state();
    }
    else
    {
        AssertMsg(false, "problem type (" + type + ") not found");
    }
}

void Transport_Problem::
get_weak_data(shared_ptr<Energy_Discretization> &energy,
              shared_ptr<Angular_Discretization> &angular,
              shared_ptr<Solid_Geometry> &solid,
              shared_ptr<Weak_Spatial_Discretization> &spatial,
              shared_ptr<Transport_Discretization> &transport,
              shared_ptr<Weak_RBF_Sweep> &sweep)
{
    // Get energy discretization
    Energy_Discretization_Parser energy_parser;
    energy = energy_parser.parse_from_xml(input_node_.get_child("energy_discretization"));

    // Get angular discretization
    Angular_Discretization_Parser angular_parser;
    angular = angular_parser.parse_from_xml(input_node_.get_child("angular_discretization"));

    // Get materials
    Material_Parser material_parser(angular,
                                    energy);
    vector<shared_ptr<Material> > materials
        = material_parser.parse_from_xml(input_node_.get_child("materials"));

    // Get boundary sources
    Boundary_Source_Parser boundary_parser(angular,
                                           energy);
    vector<shared_ptr<Boundary_Source> > boundary_sources
        = boundary_parser.parse_from_xml(input_node_.get_child("boundary_sources"));
    
    // Get solid geometry
    Constructive_Solid_Geometry_Parser solid_parser(materials,
                                                    boundary_sources);
    shared_ptr<Constructive_Solid_Geometry> constructive_solid
        = solid_parser.parse_from_xml(input_node_.get_child("solid_geometry"));
    solid = constructive_solid;
    Assert(constructive_solid->cartesian_boundaries());
    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces
        = constructive_solid->cartesian_boundary_surfaces();
    
    // Get spatial discretization
    Weak_Spatial_Discretization_Parser spatial_parser(solid,
                                                      boundary_surfaces);
    spatial = spatial_parser.get_weak_discretization(input_node_.get_child("spatial_discretization"));

    // Get transport discretization
    transport = make_shared<Transport_Discretization>(spatial,
                                                      angular,
                                                      energy);
    
    // Get sweep
    Weak_Sweep_Parser sweep_parser(spatial,
                                   angular,
                                   energy,
                                   transport);
    sweep = sweep_parser.get_weak_rbf_sweep(input_node_.get_child("transport"));
}

void Transport_Problem::
solve_steady_state()
{
    // Get preliminaries
    shared_ptr<Energy_Discretization> energy;
    shared_ptr<Angular_Discretization> angular;
    shared_ptr<Solid_Geometry> solid;
    shared_ptr<Weak_Spatial_Discretization> spatial;
    shared_ptr<Transport_Discretization> transport;
    shared_ptr<Weak_RBF_Sweep> sweep;
    get_weak_data(energy,
                  angular,
                  solid,
                  spatial,
                  transport,
                  sweep);

    // Get solver
    XML_Node solver_node = input_node_.get_child("solver");
    string type = solver_node.get_attribute<string>("type");
    Solver_Parser solver_parser(spatial,
                                angular,
                                energy,
                                transport);
    shared_ptr<Solver> solver;
    if (type == "source_iteration")
    {
        solver = solver_parser.get_source_iteration(solver_node,
                                                    sweep);
    }
    else if (type == "krylov")
    {
        solver = solver_parser.get_krylov_steady_state(solver_node,
                                                       sweep);
    }
    else
    {
        AssertMsg(false, "solver type (" + type + ") not found");
    }

    // Solve problem
    solver->solve();
    
    // Output data
    energy->output(output_node_.append_child("energy_discretization"));
    angular->output(output_node_.append_child("angular_discretization"));
    spatial->output(output_node_.append_child("spatial_discretization"));
    transport->output(output_node_.append_child("transport_discretization"));
    solid->output(output_node_.append_child("solid_geometry"));
    sweep->output(output_node_.append_child("transport"));
    solver->output(output_node_.append_child("solver"));
}

void Transport_Problem::
solve_eigenvalue()
{
    // Get preliminaries
    shared_ptr<Energy_Discretization> energy;
    shared_ptr<Angular_Discretization> angular;
    shared_ptr<Solid_Geometry> solid;
    shared_ptr<Weak_Spatial_Discretization> spatial;
    shared_ptr<Transport_Discretization> transport;
    shared_ptr<Weak_RBF_Sweep> sweep;
    get_weak_data(energy,
                  angular,
                  solid,
                  spatial,
                  transport,
                  sweep);

    // Get solver
    XML_Node solver_node = input_node_.get_child("solver");
    string type = solver_node.get_attribute<string>("type");
    Solver_Parser solver_parser(spatial,
                                angular,
                                energy,
                                transport);
    shared_ptr<Solver> solver;
    if (type == "krylov")
    {
        solver = solver_parser.get_krylov_eigenvalue(solver_node,
                                                     sweep);
    }
    else
    {
        AssertMsg(false, "solver type (" + type + ") not found");
    }

    // Solve problem
    solver->solve();
    
    // Output data
    energy->output(output_node_.append_child("energy_discretization"));
    angular->output(output_node_.append_child("angular_discretization"));
    spatial->output(output_node_.append_child("spatial_discretization"));
    transport->output(output_node_.append_child("transport_discretization"));
    solid->output(output_node_.append_child("solid_geometry"));
    sweep->output(output_node_.append_child("transport"));
    solver->output(output_node_.append_child("solver"));
}
