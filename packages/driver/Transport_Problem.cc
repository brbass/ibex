#include "Transport_Problem.hh"

#include <iostream>

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
#include "Meshless_Sweep.hh"
#include "Meshless_Sweep_Parser.hh"
#include "Solver.hh"
#include "Solver_Parser.hh"
#include "Source_Iteration.hh"
#include "Strong_Spatial_Discretization.hh"
#include "Strong_Spatial_Discretization_Parser.hh"
#include "Timer.hh"
#include "Transport_Discretization.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Parser.hh"

using namespace std;

Transport_Problem::
Transport_Problem(XML_Node input_node,
                  XML_Node output_node,
                  bool print):
    input_node_(input_node),
    output_node_(output_node),
    print_(print)
{
}

void Transport_Problem::
solve()
{
    Timer timer;
    XML_Node problem_node = input_node_.get_child("problem");
    string type = problem_node.get_attribute<string>("type");
    string discretization_method = problem_node.get_attribute<string>("discretization",
                                                                      "weak");
    
    timer.start();
    if (type == "eigenvalue")
    {
        solve_eigenvalue(discretization_method);
    }
    else if (type == "steady_state")
    {
        solve_steady_state(discretization_method);
    }
    else
    {
        AssertMsg(false, "problem type (" + type + ") not found");
    }
    timer.stop();
    times_.emplace_back(timer.time(), "total");
    output_timing();
}

void Transport_Problem::
get_weak_data(string discretization_method,
              shared_ptr<Energy_Discretization> &energy,
              shared_ptr<Angular_Discretization> &angular,
              shared_ptr<Solid_Geometry> &solid,
              shared_ptr<Weak_Spatial_Discretization> &spatial,
              shared_ptr<Transport_Discretization> &transport,
              shared_ptr<Meshless_Sweep> &sweep)
{
    Timer timer;
    
    // Get energy discretization
    print_message("Parsing energy");
    timer.start();
    Energy_Discretization_Parser energy_parser;
    energy = energy_parser.parse_from_xml(input_node_.get_child("energy_discretization"));
    timer.stop();
    times_.emplace_back(timer.time(), "energy_initialization");

    // Get angular discretization
    print_message("Parsing angular");
    timer.start();
    Angular_Discretization_Parser angular_parser;
    angular = angular_parser.parse_from_xml(input_node_.get_child("angular_discretization"));
    timer.stop();
    times_.emplace_back(timer.time(), "angular_initialization");

    // Get materials
    print_message("Parsing materials");
    timer.start();
    Material_Parser material_parser(angular,
                                    energy);
    vector<shared_ptr<Material> > materials
        = material_parser.parse_from_xml(input_node_.get_child("materials"));
    timer.stop();
    times_.emplace_back(timer.time(), "material_initialization");

    // Get boundary sources
    print_message("Parsing boundary sources");
    timer.start();
    Boundary_Source_Parser boundary_parser(angular,
                                           energy);
    vector<shared_ptr<Boundary_Source> > boundary_sources
        = boundary_parser.parse_from_xml(input_node_.get_child("boundary_sources"));
    timer.stop();
    times_.emplace_back(timer.time(), "boundary_source_initialization");
    
    // Get solid geometry
    print_message("Parsing solid geometry");
    timer.start();
    Constructive_Solid_Geometry_Parser solid_parser(materials,
                                                    boundary_sources);
    shared_ptr<Constructive_Solid_Geometry> constructive_solid
        = solid_parser.parse_from_xml(input_node_.get_child("solid_geometry"));
    solid = constructive_solid;
    Assert(constructive_solid->cartesian_boundaries());
    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces
        = constructive_solid->cartesian_boundary_surfaces();
    timer.stop();
    times_.emplace_back(timer.time(), "solid_geometry_initialization");
    
    // Get spatial discretization
    print_message("Parsing spatial discretization");
    timer.start();
    if (discretization_method == "weak")
    {
    Weak_Spatial_Discretization_Parser spatial_parser(solid,
                                                      boundary_surfaces);
    spatial = spatial_parser.get_weak_discretization(input_node_.get_child("spatial_discretization"));
    }
    else if (discretization_method == "strong")
    {
        Strong_Spatial_Discretization_Parser spatial_parser(solid,
                                                            boundary_surfaces);
        spatial = spatial_parser.get_strong_discretization(input_node_.get_child("spatial_discretization"));
    }
    else
    {
        AssertMsg(false, "discretization method (" + discretization_method + ") not found");
    }
        
    timer.stop();
    times_.emplace_back(timer.time(), "spatial_initialization");

    // Get transport discretization
    print_message("Creating transport discretization");
    timer.start();
    transport = make_shared<Transport_Discretization>(spatial,
                                                      angular,
                                                      energy);
    timer.stop();
    times_.emplace_back(timer.time(), "transport_initialization");
    
    // Get sweep
    print_message("Parsing sweep");
    timer.start();
    Meshless_Sweep_Parser sweep_parser(spatial,
                                       angular,
                                       energy,
                                       transport);
    sweep = sweep_parser.get_meshless_sweep(input_node_.get_child("transport"));
    timer.stop();
    times_.emplace_back(timer.time(), "sweep_initialization");
}

void Transport_Problem::
solve_steady_state(string discretization_method)
{
    Timer timer;
    
    // Get preliminaries
    print_message("Parsing weak data");
    timer.start();
    shared_ptr<Energy_Discretization> energy;
    shared_ptr<Angular_Discretization> angular;
    shared_ptr<Solid_Geometry> solid;
    shared_ptr<Weak_Spatial_Discretization> spatial;
    shared_ptr<Transport_Discretization> transport;
    shared_ptr<Meshless_Sweep> sweep;
    get_weak_data(discretization_method,
                  energy,
                  angular,
                  solid,
                  spatial,
                  transport,
                  sweep);
    timer.stop();
    times_.emplace_back(timer.time(), "weak_data_initialization");

    // Get solver
    print_message("Parsing solver");
    timer.start();
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
    timer.stop();
    times_.emplace_back(timer.time(), "solver_initialization");
    
    // Solve problem
    print_message("Solving problem");
    timer.start();
    solver->solve();
    timer.stop();
    times_.emplace_back(timer.time(), "solve");
    
    // Output data
    print_message("Output data");
    timer.start();
    energy->output(output_node_.append_child("energy_discretization"));
    angular->output(output_node_.append_child("angular_discretization"));
    spatial->output(output_node_.append_child("spatial_discretization"));
    transport->output(output_node_.append_child("transport_discretization"));
    solid->output(output_node_.append_child("solid_geometry"));
    sweep->output(output_node_.append_child("transport"));
    solver->output(output_node_.append_child("solver"));
    timer.stop();
    times_.emplace_back(timer.time(), "output");
}

void Transport_Problem::
solve_eigenvalue(string discretization_method)
{
    Timer timer;
    
    // Get preliminaries
    print_message("Parsing weak data");
    timer.start();
    shared_ptr<Energy_Discretization> energy;
    shared_ptr<Angular_Discretization> angular;
    shared_ptr<Solid_Geometry> solid;
    shared_ptr<Weak_Spatial_Discretization> spatial;
    shared_ptr<Transport_Discretization> transport;
    shared_ptr<Meshless_Sweep> sweep;
    get_weak_data(discretization_method,
                  energy,
                  angular,
                  solid,
                  spatial,
                  transport,
                  sweep);
    timer.stop();
    times_.emplace_back(timer.time(), "weak_data_initialization");

    // Get solver
    print_message("Parsing solver");
    timer.start();
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
    timer.stop();
    times_.emplace_back(timer.time(), "solver_initialization");
    
    // Solve problem
    print_message("Solving problem");
    timer.start();
    solver->solve();
    timer.stop();
    times_.emplace_back(timer.time(), "solve");
    
    // Output data
    print_message("Output data");
    timer.start();
    energy->output(output_node_.append_child("energy_discretization"));
    angular->output(output_node_.append_child("angular_discretization"));
    spatial->output(output_node_.append_child("spatial_discretization"));
    transport->output(output_node_.append_child("transport_discretization"));
    solid->output(output_node_.append_child("solid_geometry"));
    sweep->output(output_node_.append_child("transport"));
    solver->output(output_node_.append_child("solver"));
    timer.stop();
    times_.emplace_back(timer.time(), "output");
}

void Transport_Problem::
output_timing()
{
    XML_Node timing_node = output_node_.append_child("timing");
    for (pair<double, string> time : times_)
    {
        timing_node.set_child_value(time.first, time.second);
    }
}

void Transport_Problem::
print_message(string message) const
{
    if (print_)
    {
        cout << endl;
        cout << "Transport Problem:  ";
        cout << message << endl;
    }
}
