#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>

#include <mpi.h>

#include <Amesos.h>
#include <Epetra_Comm.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Vector.h>

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source_Parser.hh"
#include "Cartesian_Plane.hh"
#include "Check_Equality.hh"
#include "Energy_Discretization.hh"
#include "Energy_Discretization_Parser.hh"
#include "Interpolation_Solid_Geometry.hh"
#include "Material.hh"
#include "Material_Parser.hh"
#include "Random_Number_Generator.hh"
#include "Solid_Geometry.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Parser.hh"
#include "Weight_Function.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;
namespace ce = Check_Equality;

Random_Number_Generator<double> rng(-1., 1., 4109);

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
    
    return make_shared<Epetra_Map>(number_of_points, 0 /*index base*/, *comm);
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
        switch (weight->material_options().weighting)
        {
        case Weight_Function::Material_Options::Weighting::POINT:
        {
            vector<double> const v_b = weight->v_b();
            mat->InsertGlobalValues(i, number_of_basis_functions[i], &v_b[0], &basis_function_indices[0]);
            break;
        }
        case Weight_Function::Material_Options::Weighting::WEIGHT:
        {
            vector<double> const iv_b_w = weight->iv_b_w();
            mat->InsertGlobalValues(i, number_of_basis_functions[i], &iv_b_w[0], &basis_function_indices[0]);
            break;
        }
        default:
            AssertMsg(false, "weighting type not implemented");
        }
        
    }
    mat->FillComplete();
    
    return mat;
}

shared_ptr<Epetra_Vector> get_rhs(shared_ptr<Weak_Spatial_Discretization> spatial,
                                  shared_ptr<Epetra_Map> map)

{
    int number_of_points = spatial->number_of_points();
    
    shared_ptr<Epetra_Vector> vec
        = make_shared<Epetra_Vector>(*map);
    
    for (int i = 0; i < number_of_points; ++i)
    {
        int num_entries = 1;
        vector<int> global_index = {i};
        shared_ptr<Weight_Function> weight = spatial->weight(i);
        vector<double> const data = weight->material()->internal_source()->data();
        switch (weight->material_options().weighting)
        {
        case Weight_Function::Material_Options::Weighting::POINT:
        {
            (*vec)[i] = data[0] / weight->iv_w()[0];
            break;
        }
        case Weight_Function::Material_Options::Weighting::WEIGHT:
        {
            (*vec)[i] = data[0];
            break;
        }
        default:
            AssertMsg(false, "weighting option not found");
        }
    }
    
    return vec;
}

shared_ptr<Epetra_LinearProblem> get_problem(shared_ptr<Epetra_CrsMatrix> mat,
                                             shared_ptr<Epetra_Vector> lhs,
                                             shared_ptr<Epetra_Vector> rhs)
{
    return make_shared<Epetra_LinearProblem>(mat.get(),
                                             lhs.get(),
                                             rhs.get());
}

shared_ptr<Amesos_BaseSolver*> get_solver(shared_ptr<Epetra_LinearProblem> problem)
{
    Amesos factory;
    return make_shared<Amesos_BaseSolver*>(factory.Create("Klu", *problem));
}

vector<double> convert_to_vector(shared_ptr<Epetra_Vector> lhs)
{
    int number_of_points = lhs->MyLength();
    vector<double> solution(number_of_points);
    
    for (int i = 0; i < number_of_points; ++i)
    {
        solution[i] = (*lhs)[i];
    }

    return solution;
}

double get_value(shared_ptr<Weak_Spatial_Discretization> spatial,
                 vector<double> const &position,
                 vector<double> const &coefficients)
{
    return spatial->expansion_value(position,
                                    coefficients);
}

int test_interpolation(int dimension,
                       function<double(vector<double>)> const &source,
                       function<double(int, vector<double>)> const &d_source,
                       XML_Node input_node,
                       string description)
{
    int checksum = 0;
    
    shared_ptr<Weak_Spatial_Discretization> spatial
        = get_spatial(dimension,
                      source,
                      input_node);
    shared_ptr<Epetra_Comm> comm = get_comm();
    shared_ptr<Epetra_Map> map = get_map(spatial,
                                         comm);
    shared_ptr<Epetra_CrsMatrix> mat = get_matrix(spatial,
                                                  map);
    shared_ptr<Epetra_Vector> rhs = get_rhs(spatial,
                                            map);
    shared_ptr<Epetra_Vector> lhs
        = make_shared<Epetra_Vector>(*map);
    shared_ptr<Epetra_LinearProblem> problem = get_problem(mat,
                                                           lhs,
                                                           rhs);
    shared_ptr<Amesos_BaseSolver*> solver = get_solver(problem);
    (*solver)->SymbolicFactorization();
    (*solver)->NumericFactorization();
    (*solver)->Solve();

    int number_of_points = spatial->number_of_points();
    vector<double> coefficients = convert_to_vector(lhs);
    
    // Check that collocation value is equal to the original result
    {
        vector<double> values;
        spatial->collocation_values(coefficients,
                                    values);
        vector<double> expected_values(number_of_points);
        for (int i = 0; i < number_of_points; ++i)
        {
            expected_values[i] = source(spatial->weight(i)->position());
        }

        if (!ce::approx(values, expected_values, 1e-12))
        {
            checksum += 1;
            cout << "collocation failed for (" + description + ")";
        }
    }
    // Check some random points
    {
        int number_of_tests = 100;
        for (int i = 0; i < number_of_tests; ++i)
        {
            vector<double> position = rng.vector(dimension);
            
            double value = spatial->expansion_value(position,
                                                    coefficients);
            double expected_value = source(position);

            if (!ce::approx(value, expected_value, 1e-8))
            {
                checksum += 1;
                cout << "interp failed for (" + description + ") in test ";
                cout << i;
                cout << "; calculated: ";
                cout << value;
                cout << "; expected: ";
                cout << expected_value;
                cout << "; error: ";
                cout << value - expected_value;
                cout << endl;
            }
        }
    }
    
    return checksum;
}

int run_tests(string input_folder)
{
    int checksum = 0;
    
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
                = [](vector<double> const &/*position*/)
                {
                    return 1.;
                };
            function<double(int, vector<double>)> d_source
                = [](int /*dimension*/,
                     vector<double> const &/*position*/)
                {
                    return 0.;
                };
            
            checksum += test_interpolation(dimension,
                                           source,
                                           d_source,
                                           input_node,
                                           "constant");
        }
        
        // Test linear function
        {
            function<double(vector<double>)> source
                = [](vector<double> const &position)
                {
                    return 2. * position[0] + 3. * position[1];
                };
            function<double(int, vector<double>)> d_source
                = [](int dimension,
                     vector<double> const &position)
                {
                    switch (dimension)
                    {
                    case 0:
                        return 2.;
                    case 1:
                        return 3.;
                    }
                };
            
            checksum += test_interpolation(dimension,
                                           source,
                                           d_source,
                                           input_node,
                                           "linear");
        }
    }

    return checksum;
}

// Temporary function
int check_basis(string input_folder)
{
    int dimension = 2;
    string input_filename = input_folder + "/temp_test_5.xml";
    
    XML_Document input_file(input_filename);
    XML_Node input_node = input_file.get_child("input");
    function<double(vector<double>)> source
        = [](vector<double> const &/*position*/)
        {
            return 1.;
        };
    shared_ptr<Weak_Spatial_Discretization> spatial
        = get_spatial(dimension,
                      source,
                      input_node);
    
    int number_of_points = spatial->number_of_points();
    for (int i = 0; i < number_of_points; ++i)
    {
        shared_ptr<Weight_Function> weight = spatial->weight(i);
        
        double const iv_w = weight->iv_w()[0];
        vector<double> const iv_b_w = weight->iv_b_w();
        
        double sum = 0.;
        for (double val : iv_b_w)
        {
            sum += val;
        }
        
        int w = 16;
        cout << setw(w) << iv_w;
        cout << setw(w) << sum;
        cout << setw(w) << iv_w - sum;
        cout << endl;
    }
}

int main(int argc, char **argv)
{
    int checksum = 0;

    MPI_Init(&argc, &argv);
    
    if (argc != 2)
    {
        cerr << "usage: tst_Interpolation [input_folder]" << endl;
        return 1;
    }
    
    string input_folder = argv[1];
    input_folder += "/";

    // check_basis(input_folder);
    run_tests(input_folder);
    
    MPI_Finalize();
    
    return checksum;
}
