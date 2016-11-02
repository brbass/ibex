#include <iomanip>
#include <iostream>
#include <memory>
#include <functional>

#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source_Parser.hh"
#include "Check_Equality.hh"
#include "Energy_Discretization_Parser.hh"
#include "Material_Parser.hh"
#include "RBF_Discretization.hh"
#include "Spatial_Discretization_Parser.hh"
#include "Trilinos_Dense_Solve.hh"
#include "Vector_Functions.hh"

using namespace std;
namespace ce = Check_Equality;

shared_ptr<RBF_Discretization> get_rbf_discretization(string xml_input_filename)
{
    pugi::xml_document input_document;

    if (!input_document.load_file(xml_input_filename.c_str()))
    {
        cout << "Could not open xml input file \"" + xml_input_filename + "\"" << endl;
        return shared_ptr<RBF_Discretization>();
    }

    pugi::xml_node input_file = input_document.child("input");

    Angular_Discretization_Parser angular_parser(input_file);
    shared_ptr<Angular_Discretization> angular = angular_parser.get_ptr();
    
    Energy_Discretization_Parser energy_parser(input_file);
    shared_ptr<Energy_Discretization> energy = energy_parser.get_ptr();
    
    Material_Parser material_parser(input_file,
                                    angular,
                                    energy);
    vector<shared_ptr<Material> > materials = material_parser.get_ptr();
    
    Boundary_Source_Parser boundary_parser(input_file,
                                           angular,
                                           energy);
    vector<shared_ptr<Boundary_Source> > boundary_sources = boundary_parser.get_ptr();

    Spatial_Discretization_Parser spatial_parser(input_file,
                                                 materials,
                                                 boundary_sources);
    shared_ptr<Spatial_Discretization> spatial = spatial_parser.get_ptr();
    
    return dynamic_pointer_cast<RBF_Discretization>(spatial);
}

double sin_function(vector<double> const &x)
{
    return sin(M_PI * x[0]) * sin(M_PI * x[1]);
}

vector<double> gradient_sin_function(vector<double> const &x)
{
    return {M_PI * cos(M_PI * x[0]) * sin(M_PI * x[1]), 
            M_PI * sin(M_PI * x[0]) * cos(M_PI * x[1])};
}

double flat_function(vector<double> const &x)
{
    return 1.;
}

vector<double> gradient_flat_function(vector<double> const &x)
{
    return {0., 0.};
}

double linear_function(vector<double> const &x)
{
    return x[0] * x[1];
}

vector<double> gradient_linear_function(vector<double> const &x)
{
    return {x[1], x[0]};
}

vector<double> get_basis_matrix(shared_ptr<RBF_Discretization> rbf)
{
    int number_of_points = rbf->number_of_points();
    
    vector<double> matrix_data(number_of_points * number_of_points);

    for (int i = 0; i < number_of_points; ++i)
    {
        for (int j = 0; j < number_of_points; ++j)
        {
            matrix_data[j + number_of_points * i] = rbf->basis(i, j, 0, 0);
        }
    }
    
    return matrix_data;
}

int test_interpolation(string xml_input_filename)
{
    int checksum = 0;

    double tolerance = 1e-4;
    double gradient_tolerance = 1e-1;

    cout << "Running input file " << xml_input_filename << endl;
    
    // Get RBF Discretization
    shared_ptr<RBF_Discretization> rbf = get_rbf_discretization(xml_input_filename);
    if (!rbf)
    {
        checksum += 1;
        return checksum;
    }
    int dimension = rbf->dimension();
    int number_of_points = rbf->number_of_points();

    // Get matrix data
    vector<double> const matrix_data = get_basis_matrix(rbf);
    
    Trilinos_Dense_Solve solver;

    // List of functions to use as RHS
    int num_funcs = 3;
    vector<function<double(vector<double>)> > funcs
        = {flat_function, linear_function, sin_function};
    vector<function<vector<double>(vector<double>)> > gradient_funcs
        = {gradient_flat_function, gradient_linear_function, gradient_sin_function};
    vector<string> function_descriptions = {"f(x, y) = 1",
                                            "f(x, y) = xy",
                                            "f(x, y) = sin(x) sin(x)"};
    
    for (int f = 0; f < num_funcs; ++f)
    {
        cout << "    Function type:\t" << function_descriptions[f] << endl;
        function<double(vector<double>)> func = funcs[f];
        function<vector<double>(vector<double>)> gradient_func = gradient_funcs[f];
        
        // Find coefficients for each RHS function
        vector<double> rhs_data(number_of_points);
        vector<double> gradient_data(number_of_points * dimension);
        for (int i = 0; i < number_of_points; ++i)
        {
            rhs_data[i] = func(rbf->rbf_point(i)->position());

            vector<double> gradient = gradient_func(rbf->rbf_point(i)->position());
            for (int d = 0; d < dimension; ++d)
            {
                gradient_data[d + dimension * i] = gradient[d];
            }
        }
        
        // Solve system of equations
        vector<double> solution_data(number_of_points);
        vector<double> matrix_temp = matrix_data;
        vector<double> rhs_temp = rhs_data;

        solver.epetra_solve(matrix_temp,
                            rhs_temp,
                            solution_data,
                            number_of_points);
        
        // Get list of values and derivative values
        vector<double> values(number_of_points, 0);
        vector<double> gradient_values(number_of_points * dimension, 0);
        
        for (int i = 0; i < number_of_points; ++i)
        {
            for (int j = 0; j < number_of_points; ++j)
            {
                values[i] += solution_data[j] * rbf->basis(i, j, 0, 0);
                
                vector<double> gradient = rbf->gradient_basis(i, j, 0, 0);
                for (int d = 0; d < dimension; ++d)
                {
                    gradient_values[d + dimension * i] += solution_data[j] * gradient[d];
                }
            }
        }
        
        bool printed_header = false;
        for (int i = 0; i < number_of_points; ++i)
        {
            if (!ce::approx(values[i], rhs_data[i], tolerance))
            {
                int w = 16;
                if (!printed_header)
                {
                    cout << "\tValues incorrect:" << endl;
                    cout << "\t";
                    cout << setw(w) << "point";
                    cout << setw(w) << "calculated";
                    cout << setw(w) << "expected";
                    cout << endl;
                    printed_header = true;
                }
                cout << "\t";
                cout << setw(w) << i;
                cout << setw(w) << values[i];
                cout << setw(w) << rhs_data[i];
                cout << endl;
                
                checksum += 1;
            }
        }
        
        printed_header = false;
        for (int i = 0; i < number_of_points; ++i)
        {
            for (int d = 0; d < dimension; ++d)
            {
                if (!ce::approx(gradient_values[d + dimension * i], gradient_data[d + dimension * i], gradient_tolerance))
                {
                    int w = 16;
                    if (!printed_header)
                    {
                        cout << "\tGradient values incorrect:" << endl;
                        cout << "\t";
                        cout << setw(w) << "point";
                        cout << setw(w) << "dimension";
                        cout << setw(w) << "calculated";
                        cout << setw(w) << "expected";
                        cout << endl;
                        printed_header = true;
                    }
                    cout << "\t";
                    cout << setw(w) << i;
                    cout << setw(w) << d;
                    cout << setw(w) << gradient_values[d + dimension * i];
                    cout << setw(w) << gradient_data[d + dimension * i];
                    cout << endl;
                    
                    checksum += 1;
                }
            }
        }
    }
    cout << endl;

    return checksum;
}

int main(int argc, char **argv)
{
    int checksum = 0;

    if (argc != 2)
    {
        cerr << "usage: tst_RBF_Interpolation [input_folder]" << endl;
        return 1;
    }
    
    string input_folder = argv[1];
    input_folder += "/";
    vector<string> input_files = {"two_region_interpolation_cartesian.xml",
                                  "two_region_interpolation_optical.xml"};
    
    for (string const &input_file : input_files)
    {
        checksum += test_interpolation(input_folder + input_file);
    }

    return checksum;
}
