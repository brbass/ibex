#include <functional>

#include <mpi.h>

#include <Amesos.h>
#include <Epetra_Comm.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_MultiVector.h>

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source_Parser.hh"
#include "Cartesian_Plane.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Interpolation_Solid_Geometry.hh"
#include "Material.hh"
#include "Material_Parser.hh"
#include "Solid_Geometry.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "Weight_Function.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;

shared_ptr<Weak_Spatial_Discretization> get_spatial(int dimension,
                                                    function<double(vector<double>)> const &source,
                                                    XML_Node input_node)
{
    // Get angular discretization
    Angular_Discretization_Parser angular_parser;
    shared_ptr<Angular_Discretization> angular
        = angular_parser.parse_from_xml(input_node.get_child("angular_discretization"));

    // Get energy discretization
    Energy_Discretization_Parser energy_parser;
    shared_ptr<Energy_Discretization> energy
        = energy_parser.parse_from_xml(input_node.get_child("energy_discretization"));

    // Get boundary sources
    Boundary_Source_Parser boundary_parser(angular,
                                           energy);
    vector<shared_ptr<Boundary_Source> > boundary_sources
        = boundary_parser.parse_from_xml(input_node.get_child("boundary_sources"));

    // Get solid geometry
    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces(2 * dimension);
    for (int d = 0; d < dimension; ++d)
    {
        for (int i = 0; i < 2; ++i)
        {
            int index = d + dimension * i;
            double position = i == 0 ? -1 : 1;
            double normal = i == 0 ? -1 : 1;
            boundary_surfaces[index]
                = make_shared<Cartesian_Plane>(index,
                                               dimension,
                                               Surface::Surface_Type::BOUNDARY,
                                               d,
                                               position,
                                               normal);
            boundary_surfaces[index]->set_boundary_source(boundary_sources[0]);
        }
                                               
    }
    shared_ptr<Solid_Geometry> solid_geometry
        = make_shared<Analytic_Solid_Geometry>(dimension,
                                               angular,
                                               energy,
                                               source);
    
    // Parser for basis and weight functions
    Weak_Spatial_Discretization_Parser spatial_parser(solid_geometry,
                                                      boundary_surfaces);
    return spatial_parser.get_weak_discretization(input_node.get_child("spatial_discretization"));
}

shared_ptr<Epetra_Comm> get_comm()
{
    return make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
}

shared_ptr<Epetra_Map> get_map(shared_ptr<Weak_Spatial_Discretization> spatial,
                               shared_ptr<Epetra_Comm> comm)
{
    int number_of_points = spatial->number_of_points();
    
    return make_shared<Epetra_Map>(number_of_points, index_base_, *comm);
}

shared_ptr<Epetra_CrsMatrix> get_matrix(shared_ptr<Weak_Spatial_Discretization> spatial,
                                        shared_ptr<Epetra_Map> map)
{
    int number_of_points = spatial->number_of_points();
    vector<int> const number_of_basis_functions = spatial->number_of_basis_functions();
    shared_ptr<Epetra_CrsMatrix> mat
        = make_shared<Epetra_CrsMatrix>(Copy, *map, &number_of_basis_functions[0], true);
    
    for (int i = 0; i < number_of_points; ++i)
    {
        shared_ptr<Weight_Function> weight = spatial->weight(i);
        vector<int> const basis_function_indices = weight->basis_function_indices();
        vector<double> vals(number_of_basis_functions[i]);
        switch(weight->material_options().weighting)
        {
        case Weight_Function::Material_Options::Weighting::POINT:
        {
            vector<double> const v_b = weight->v_b();
            mat->InsertGlobalValues(i, number_of_basis_functions[i], &v_b[0], &basis_function_indices[0]);
        }
        case Weight_Function::Material_Options::Weighting::WEIGHT:
        {
            vector<double> const iv_b_w = weight->iv_b_w();
            mat->InsertGlobalValues(i, number_of_basis_functions[i], &iv_b_w[0], &basis_function_indices[0]);
        }
        default:
            AssertMsg(false, "weighting type not implemented");
        }
        
    }
    mat->FillComplete();
    
    return mat;
}

shared_ptr<Epetra_Vector> get_rhs(shared_ptr<Map> map,
                                  shared_ptr<Weak_Spatial_Discretization> spatial)
{
    int number_of_points = spatial->number_of_points();

    shared_ptr<Epetra_Vector> vec
        = make_shared<Epetra_Vector>(*map);
    
    for (int i = 0; i < number_of_points; ++i)
    {
        int num_entries = 1;
        vector<int> global_index = {i};
        vector<double> const data = spatial->weight(i)->material()->internal_source()->data();
        vec->ReplaceGlobalValues(num_entries,
                                 &data[0],
                                 &global_index[0]);
    }
}

shared_ptr<Epetra_LinearProblem> get_problem(shared_ptr<Epetra_CrsMatrix> mat,
                                             shared_ptr<Epetra_Vector> lhs,
                                             shared_ptr<Epetra_Vector> rhs)
{
    return make_shared<Epetra_LinearProblem>(matrix_.get(),
                                             lhs_.get(),
                                             rhs_.get());
}

shared_ptr<Amesos_BaseSolver> get_solver(shared_ptr<Epetra_LinearProblem> problem)
{
    return make_shared<Amesos_BaseSolver>(factory_.Create("Klu", *problem));
}

int test_interpolation(int dimension,
                       function<double(vector<double>)> const &source,
                       function<double(vector<double>)> const &d_source,
                       XML_Node input_node)
{
    int checksum = 0;
    
    shared_ptr<Weak_Spatial_Discretization> spatial
        = get_spatial(dimension,
                      source,
                      input_node);

    
    
    
    
    return checksum;
}

int main(int argc, char **argv)
{
    int checksum = 0;

    if (argc != 2)
    {
        cerr << "usage: tst_Interpolation [input_folder]" << endl;
        return 1;
    }
    
    string input_folder = argv[1];
    input_folder += "/";
    vector<string> input_filenames
        = {input_folder + "/mls_interpolation.xml"};
    
    for (string input_filename : input_filenames)
    {
        // Get XML document
        XML_Document input_file(input_filename);
        XML_Node input_node = input_file.get_child("input");
        int dimension = input_node.get_child("angular_discretization").get_child_value<int>("dimension");
        
        // Test constant
        {
            function<double(vector<double>)> source
                = [](vector<double> const &position)
                {
                    return 1.;
                };
            function<double(vector<double>)> d_source
                = [](vector<double> const &position)
                {
                    return 0.;
                };
            
            checksum += test_interpolation(dimension,
                                           source,
                                           d_source,
                                           input_node);
        }
    }
    
    return checksum;
}