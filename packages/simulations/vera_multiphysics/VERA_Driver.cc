#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source_Parser.hh"
#include "Cartesian_Plane.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Krylov_Eigenvalue.hh"
#include "LDFE_Quadrature.hh"
#include "Material_Parser.hh"
#include "Solver.hh"
#include "Solver_Parser.hh"
#include "Transport_Discretization.hh"
#include "VERA_Solid_Geometry.hh"
#include "Weak_RBF_Sweep.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "Weak_Sweep_Parser.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;

void test_problem(string filename)
{
    XML_Document input_file(filename);
    XML_Node input_node = input_file.get_child("input");
    
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
        = make_shared<VERA_Solid_Geometry>(false,
                                           [](vector<double> const &){return 600;},
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

    string filename = argv[1];

    test_problem(filename);
        
    // Close MPI
    MPI_Finalize();

    return 0;
}
