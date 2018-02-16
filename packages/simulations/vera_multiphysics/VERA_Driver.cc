#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source_Parser.hh"
#include "Cartesian_Plane.hh"
#include "Cross_Section.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Heat_Transfer_Integration.hh"
#include "Heat_Transfer_Factory.hh"
#include "Heat_Transfer_Solve.hh"
#include "Heat_Transfer_Solution.hh"
#include "Krylov_Eigenvalue.hh"
#include "LDFE_Quadrature.hh"
#include "Material.hh"
#include "Material_Parser.hh"
#include "Quadrature_Rule.hh"
#include "Solver.hh"
#include "Solver_Parser.hh"
#include "Transport_Discretization.hh"
#include "VERA_Heat_Data.hh"
#include "VERA_Solid_Geometry.hh"
#include "VERA_Transport_Result.hh"
#include "Weak_RBF_Sweep.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "Weak_Sweep_Parser.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;

shared_ptr<VERA_Transport_Result>
run_transport(XML_Node input_node,
              shared_ptr<VERA_Temperature> temperature)
{
    // Get energy discretization
    Energy_Discretization_Parser energy_parser;
    shared_ptr<Energy_Discretization> energy
        = energy_parser.parse_from_xml(input_node.get_child("energy_discretization"));

    // Get angular discretization
    Angular_Discretization_Parser angular_parser;
    shared_ptr<Angular_Discretization> angular
        = angular_parser.parse_from_xml(input_node.get_child("angular_discretization"));

    // Get materials
    Material_Parser material_parser(angular,
                                    energy);
    bool include_ifba = input_node.get_child("materials").get_attribute<bool>("include_ifba");
    vector<shared_ptr<Material> > materials
        = material_parser.parse_from_xml(input_node.get_child("materials"));
    
    // Get boundary source
    Boundary_Source_Parser boundary_parser(angular,
                                           energy);
    vector<shared_ptr<Boundary_Source> > boundary_sources
        = boundary_parser.parse_from_xml(input_node.get_child("boundary_sources"));
    Assert(boundary_sources.size() == 1);
    
    // Get solid geometry
    shared_ptr<VERA_Solid_Geometry> solid
        = make_shared<VERA_Solid_Geometry>(include_ifba,
                                           temperature,
                                           angular,
                                           energy,
                                           materials,
                                           boundary_sources[0]);
    
    // Get boundary surfaces
    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces
        = solid->cartesian_boundary_surfaces();
    
    // Get spatial discretization
    Weak_Spatial_Discretization_Parser spatial_parser(solid,
                                                      boundary_surfaces);
    shared_ptr<Weak_Spatial_Discretization> spatial
        = spatial_parser.get_weak_discretization(input_node.get_child("spatial_discretization"));
    
    // Get transport discretization
    shared_ptr<Transport_Discretization> transport
        = make_shared<Transport_Discretization>(spatial,
                                                angular,
                                                energy);

    // Get sweep
    Weak_Sweep_Parser sweep_parser(spatial,
                                   angular,
                                   energy,
                                   transport);
    shared_ptr<Weak_RBF_Sweep> sweep
        = sweep_parser.get_weak_rbf_sweep(input_node.get_child("transport"));
    
    // Get solver
    Solver_Parser solver_parser(spatial,
                                angular,
                                energy,
                                transport);
    shared_ptr<Solver> solver
        = solver_parser.get_krylov_eigenvalue(input_node.get_child("solver"),
                                              sweep);
    solver->solve();

    cout << "eigenvalue: " << "\t" << solver->result()->k_eigenvalue << endl;
    
    return make_shared<VERA_Transport_Result>(solid,
                                              angular,
                                              energy,
                                              spatial,
                                              solver->result());
}

shared_ptr<VERA_Temperature> 
run_heat(XML_Node input_node,
         shared_ptr<VERA_Transport_Result> result)
{
    // Get solid geometry
    double length = 0.475;
    shared_ptr<Solid_Geometry> solid;
    vector<shared_ptr<Cartesian_Plane> > surfaces;
    Heat_Transfer_Factory factory;
    factory.get_solid_1d(length,
                         solid,
                         surfaces);

    // Get weak spatial discretization
    Weak_Spatial_Discretization_Parser spatial_parser(solid,
                                                      surfaces);
    shared_ptr<Weak_Spatial_Discretization> spatial
        = spatial_parser.get_weak_discretization(input_node.get_child("spatial_discretization"));
    
    // Get heat transfer data
    shared_ptr<VERA_Heat_Data> data
        = make_shared<VERA_Heat_Data>(result);

    // Get heat transfer integration
    shared_ptr<Heat_Transfer_Integration_Options> integration_options
        = make_shared<Heat_Transfer_Integration_Options>();
    integration_options->geometry = Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_1D;
    shared_ptr<Heat_Transfer_Integration> integration
        = make_shared<Heat_Transfer_Integration>(integration_options,
                                                 data,
                                                 spatial);
    
    // Get heat transfer solver
    shared_ptr<Heat_Transfer_Solve> solver
        = make_shared<Heat_Transfer_Solve>(integration,
                                           spatial);
    shared_ptr<Heat_Transfer_Solution> solution
        = solver->solve();
    
    return make_shared<VERA_Temperature>([solution](vector<double> const & position) -> double
                                         {
                                             double radius = sqrt(position[0] * position[0] + position[1] * position[1]);
                                             if (radius > 0.475)
                                             {
                                                 return 600.0;
                                             }
                                             else
                                             {
                                                 return solution->solution({radius});
                                             }
                                         });
}

void run_test(XML_Node input_node)
{
    // Get initial temperature
    shared_ptr<VERA_Temperature> temperature
        = make_shared<VERA_Temperature>([](vector<double> const &){return 120;});

    // Run transport calculation
    shared_ptr<VERA_Transport_Result> result;

    for (int i = 0; i < 4; ++i)
    {
        result
            = run_transport(input_node.get_child("transport"),
                            temperature);
        
        // Run heat transfer calculation
        temperature
            = run_heat(input_node.get_child("heat"),
                       result);
        
        for (double r = 0; r < 0.476; r += 0.025)
        {
            vector<double> position = {r, 0};
            
            cout << r << "\t" << (*temperature)(position) << endl;
        }
    }
}

int main(int argc, char **argv)
{
    // Initialize MPI
    MPI_Init(&argc, &argv);

    if (argc != 2)
    {
        cerr << "need input file" << endl;
        return 1;
    }

    // Get base XML node
    string filename = argv[1];
    XML_Document input_file(filename);
    XML_Node input_node = input_file.get_child("input");

    run_test(input_node);
    
    // Close MPI
    MPI_Finalize();

    return 0;
}
