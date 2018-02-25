#include "Manufactured_Problem.hh"

#include <iostream>

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Parser.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Constructive_Solid_Geometry_Parser.hh"
#include "Conversion.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Integration_Mesh.hh"
#include "Krylov_Eigenvalue.hh"
#include "Krylov_Steady_State.hh"
#include "Manufactured_Integral_Operator.hh"
#include "Manufactured_Parser.hh"
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

Manufactured_Problem::
Manufactured_Problem(XML_Node input_node,
                     XML_Node output_node,
                     bool print):
    input_node_(input_node),
    output_node_(output_node),
    print_(print)
{
}

void Manufactured_Problem::
solve()
{
    Timer timer;
    XML_Node problem_node = input_node_.get_child("problem");
    string discretization_method = problem_node.get_attribute<string>("discretization",
                                                                      "weak");
    timer.start();
    solve_steady_state(discretization_method);
    timer.stop();
    times_.emplace_back(timer.time(), "total");
    output_timing();
}

void Manufactured_Problem::
get_weak_data(string discretization_method,
              shared_ptr<Energy_Discretization> &energy,
              shared_ptr<Angular_Discretization> &angular,
              shared_ptr<Manufactured_Solution> &solution,
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

    // Get solution and geometry
    print_message("Parsing solution and solid geometry");
    timer.start();
    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces;
    Manufactured_Parser manufactured_parser(angular,
                                            energy);
    manufactured_parser.parse_from_xml(input_node_.get_child("manufactured"),
                                       solution,
                                       solid,
                                       boundary_surfaces);
    timer.stop();
    times_.emplace_back(timer.time(), "solution_initialization");

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

void Manufactured_Problem::
solve_steady_state(string discretization_method)
{
    Timer timer;
    
    // Get preliminaries
    print_message("Parsing weak data");
    timer.start();
    shared_ptr<Energy_Discretization> energy;
    shared_ptr<Angular_Discretization> angular;
    shared_ptr<Manufactured_Solution> solution;
    shared_ptr<Solid_Geometry> solid;
    shared_ptr<Weak_Spatial_Discretization> spatial;
    shared_ptr<Transport_Discretization> transport;
    shared_ptr<Meshless_Sweep> sweep;
    get_weak_data(discretization_method,
                  energy,
                  angular,
                  solution,
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
    
    // Compare results
    vector<string> error_types;
    vector<vector<double> > errors;
    XML_Node manufactured_node = input_node_.get_child("manufactured");
    for (XML_Node error_node = manufactured_node.get_child("error");
         error_node;
         error_node = error_node.get_sibling("error",
                                             false))
    {
        vector<double> error;
        
        // Get integration options
        int dimension = spatial->dimension();
        std::shared_ptr<Integration_Mesh_Options> integration_options
            = make_shared<Integration_Mesh_Options>();
        integration_options->initialize_from_weak_options(spatial->options());
        integration_options->adaptive_quadrature
            = error_node.get_attribute<bool>("adaptive_quadrature",
                                             false);
        if (integration_options->adaptive_quadrature)
        {
            integration_options->minimum_radius_ordinates
                = error_node.get_attribute<bool>("minimum_radius_ordinates");
        }
        
        integration_options->integration_ordinates
            = error_node.get_attribute<int>("integration_ordinates",
                                            integration_options->integration_ordinates);
        integration_options->limits
            = error_node.get_child_matrix<double>("limits",
                                                  dimension,
                                                  2,
                                                  integration_options->limits);
        integration_options->dimensional_cells
            = error_node.get_child_vector<int>("dimensional_cells",
                                               dimension,
                                               integration_options->dimensional_cells);

        // Get error operator
        Manufactured_Integral_Operator::Options options;
        string norm = error_node.get_attribute<string>("norm");
        options.norm = options.norm_conversion()->convert(norm);
        shared_ptr<Manufactured_Integral_Operator> oper
            = make_shared<Manufactured_Integral_Operator>(options,
                                                          integration_options,
                                                          spatial,
                                                          angular,
                                                          energy,
                                                          solution);

        // Find error
        error = solver->result()->coefficients;

        (*oper)(error);

        // Store results
        error_types.push_back(norm);
        errors.push_back(error);
    }
    
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
    if (errors.size() > 0)
    {
        XML_Node error_node = output_node_.append_child("error");
        for (int i = 0; i < errors.size(); ++i)
        {
            error_node.set_child_vector<double>(errors[i], error_types[i], "group-moment-cell");
        }
    }
    timer.stop();
    times_.emplace_back(timer.time(), "output");
}

void Manufactured_Problem::
output_timing()
{
    XML_Node timing_node = output_node_.append_child("timing");
    for (pair<double, string> time : times_)
    {
        timing_node.set_child_value(time.first, time.second);
    }
}

void Manufactured_Problem::
print_message(string message) const
{
    if (print_)
    {
        cout << endl;
        cout << "Manufactured Problem:  ";
        cout << message << endl;
    }
}
