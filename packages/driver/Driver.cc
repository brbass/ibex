#include "Driver.hh"

#include "pugixml.hh"

#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source_Parser.hh"
#include "Check.hh"
#include "Energy_Discretization_Parser.hh"
#include "Material_Parser.hh"
#include "Solver_Parser.hh"
#include "Spatial_Discretization_Parser.hh"
#include "Sweep_Parser.hh"
#include "Timer.hh"
#include "Transport_Discretization.hh"
#include "Transport_Problem_Parser.hh"
#include "Vector_Operator.hh"

using namespace std;

Driver::
Driver(string xml_in):
    xml_in_(xml_in)
{
    xml_out_ = xml_in_ + ".out";
    
    run_problem();
}

void Driver::
run_problem()
{
    // folder_ = filename.substr(0, filename.find_last_of("/\\") + 1);
    
    pugi::xml_document input_document;
    
    if (!input_document.load_file(xml_in_.c_str()))
    {
        AssertMsg(false, "Could not open xml input file \"" + xml_in_ + "\"");
    }
    
    pugi::xml_node input_file = input_document.child("input");

    Timer total_timer;
    total_timer.start();

    Timer timer;
    vector<double> time(0);
    vector<string> time_desc(0);
    
    // Angular and energy discretization

    timer.start();
    
    shared_ptr<Angular_Discretization_Parser> angular_parser
        = make_shared<Angular_Discretization_Parser>(input_file);
    shared_ptr<Angular_Discretization> angular = angular_parser->get_ptr();

    timer.stop();
    add_time(timer.time(),
             "angular_parser");
    
    timer.start();
    
    shared_ptr<Energy_Discretization_Parser> energy_parser
        = make_shared<Energy_Discretization_Parser>(input_file);
    shared_ptr<Energy_Discretization> energy = energy_parser->get_ptr();

    timer.stop();
    add_time(timer.time(),
             "energy_parser");
    
    // Physical data

    timer.start();
    
    shared_ptr<Boundary_Source_Parser> boundary_parser
        = make_shared<Boundary_Source_Parser>(input_file,
                                              angular,
                                              energy);
    vector<shared_ptr<Boundary_Source> > boundary_sources = boundary_parser->get_ptr();

    timer.stop();
    add_time(timer.time(),
             "boundary_source_parser");
    
    timer.start();
    
    shared_ptr<Material_Parser> material_parser
        = make_shared<Material_Parser>(input_file,
                                       angular,
                                       energy);
    vector<shared_ptr<Material> > materials = material_parser->get_ptr();

    timer.stop();
    add_time(timer.time(),
             "material_parser");
    
    // Spatial Discretization

    timer.start();
    
    shared_ptr<Spatial_Discretization_Parser> spatial_parser
        = make_shared<Spatial_Discretization_Parser>(input_file,
                                                     materials,
                                                     boundary_sources);
    shared_ptr<Spatial_Discretization> spatial = spatial_parser->get_ptr();
    
    timer.stop();
    add_time(timer.time(),
             "spatial_discretization_parser");
    
    // Transport Discretization

    timer.start();
    
    shared_ptr<Transport_Discretization> transport
        = make_shared<Transport_Discretization>(spatial,
                                                angular,
                                                energy);
    timer.stop();
    add_time(timer.time(),
             "transport_discretization_parser");
    
    // Sweep operator

    timer.start();
    
    shared_ptr<Sweep_Parser> sweep_parser
        = make_shared<Sweep_Parser>(input_file,
                                    spatial,
                                    angular,
                                    energy,
                                    transport);
    shared_ptr<Sweep_Operator> sweeper = sweep_parser->get_ptr();

    timer.stop();
    add_time(timer.time(),
             "sweep_operator_parser");

    // Transport operator and solver

    timer.start();
    
    shared_ptr<Solver_Parser> solver_parser
        = make_shared<Solver_Parser>(input_file,
                                     spatial,
                                     angular,
                                     energy,
                                     transport,
                                     sweeper);
    shared_ptr<Solver> solver = solver_parser->get_ptr();

    timer.stop();
    add_time(timer.time(),
             "solver_parser");

    // Transport problem

    timer.start();
    
    shared_ptr<Transport_Problem_Parser> problem_parser
        = make_shared<Transport_Problem_Parser>(input_file,
                                                solver);
    
    shared_ptr<Transport_Problem> problem = problem_parser->get_ptr();
    
    timer.stop();
    add_time(timer.time(),
             "transport_problem_parser");

    // Solve problem

    timer.start();
    
    problem->solve();

    timer.stop();
    add_time(timer.time(),
             "solution");
    
    // Output data
    
    timer.start();
    
    pugi::xml_document output_document;
    pugi::xml_node output_file = output_document.append_child("output");

    angular->output(output_file);
    energy->output(output_file);
    spatial->output(output_file);
    sweeper->output(output_file);
    solver->output(output_file);
    problem->output(output_file);
    
    timer.stop();
    add_time(timer.time(),
             "output");
    
    total_timer.stop();
    add_time(total_timer.time(),
             "total");
    
    // Output timing
    
    pugi::xml_node timing = output_file.append_child("timing");

    for (int i = 0; i < times_.size(); ++i)
    {
        XML_Functions::append_child(timing, times_[i], times_description_[i]);
    }

    // Save output document
    
    output_document.save_file(xml_out_.c_str());
}
