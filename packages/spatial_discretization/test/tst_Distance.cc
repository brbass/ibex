#include <iostream>
#include <memory>

#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source_Parser.hh"
#include "Cartesian_Distance.hh"
#include "Check_Equality.hh"
#include "Distance.hh"
#include "Energy_Discretization_Parser.hh"
#include "Material_Parser.hh"
#include "Optical_Distance.hh"
#include "Constructive_Solid_Geometry_Parser.hh"
#include "Vector_Functions.hh"

using namespace std;

namespace ce = Check_Equality;
namespace vf = Vector_Functions;

// Only checks 2D

shared_ptr<Constructive_Solid_Geometry> get_solid_geometry(string xml_input_filename)
{
    pugi::xml_document input_document;
    
    if (!input_document.load_file(xml_input_filename.c_str()))
    {
        cout << "Could not open xml input file \"" + xml_input_filename + "\"" << endl;
        return shared_ptr<Constructive_Solid_Geometry>();
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
    
    Constructive_Solid_Geometry_Parser solid_parser(input_file,
                                       materials,
                                       boundary_sources);
    return solid_parser.get_ptr();
}

int test_distance(int test_case,
                  shared_ptr<Distance> distance_metric,
                  double const expected_distance,
                  vector<double> const &expected_gradient_distance,
                  vector<double> const &expected_double_gradient_distance,
                  vector<double> const &initial_position,
                  vector<double> const &final_position)
{
    int checksum = 0;
    int number_of_groups = 1;
    int group = 0;
    int dimension = 2;
    double tol = 1e-15;
    string description = distance_metric->description();
    
    // Check distance and mean_distance
    
    double distance = distance_metric->distance(group,
                                                final_position,
                                                initial_position);
    
    if (!ce::approx(distance, expected_distance, tol))
    {
        cout << description;
        cout << " case ";
        cout << test_case;
        cout << ": distance incorrect" << endl;
        cout << "\tactual: ";
        cout << expected_distance;
        cout << "\tcalculated: ";
        cout << distance << endl;
        cout << distance - expected_distance;
        checksum += 1;
    }
    
    // Check d_distance

    for (int d = 0; d < dimension; ++d)
    {
        double d_distance = distance_metric->d_distance(group,
                                                        d,
                                                        final_position,
                                                        initial_position);
        
        if (!ce::approx(d_distance, expected_gradient_distance[d], tol))
        {
            cout << description;
            cout << " case ";
            cout << test_case;
            cout << ": d_distance incorrect" << endl;
            checksum += 1;
        }
    }

    // Check dd_distance

    for (int d = 0; d < dimension; ++d)
    {
        double dd_distance = distance_metric->dd_distance(group,
                                                          d,
                                                          final_position,
                                                          initial_position);
        
        if (!ce::approx(dd_distance, expected_double_gradient_distance[d + dimension * d], tol))
        {
            cout << description;
            cout << " case ";
            cout << test_case;
            cout << ": dd_distance incorrect" << endl;
            checksum += 1;
        }
    }

    // Check gradient_distance

    vector<double> gradient_distance = distance_metric->gradient_distance(group,
                                                                          final_position,
                                                                          initial_position);
    
    if (!ce::approx(gradient_distance, expected_gradient_distance, tol))
    {
        cout << description;
        cout << " case ";
        cout << test_case;
        cout << ": gradient_distance incorrect" << endl;
        checksum += 1;
    }
    
    // Check double gradient distance

    vector<double> double_gradient_distance = distance_metric->double_gradient_distance(group,
                                                                                        final_position,
                                                                                        initial_position);

    if (!ce::approx(double_gradient_distance, expected_double_gradient_distance, tol))
    {
        cout << description;
        cout << " case ";
        cout << test_case;
        cout << ": double_gradient_distance incorrect" << endl;
        checksum += 1;
    }

    // Check laplacian distance
    
    double laplacian_distance = distance_metric->laplacian_distance(group,
                                                                    final_position,
                                                                    initial_position);

    double expected_laplacian_distance = 0;

    for (int d = 0; d < dimension; ++d)
    {
        expected_laplacian_distance += expected_double_gradient_distance[d + dimension * d];
    }
    
    if (!ce::approx(laplacian_distance, expected_laplacian_distance, tol))
    {
        cout << description;
        cout << " case ";
        cout << test_case;
        cout << ": laplacian_distance incorrect" << endl;
        checksum += 1;
    }
    
    return checksum;
}

int main(int argc, char **argv)
{
    int checksum = 0;
    
    if (argc != 2)
    {
        cerr << "usage: tst_Distance [input_folder]" << endl;
        return 1;
    }

    int const dimension = 2;
    
    string input_folder = argv[1];
    string xml_input_filepath = input_folder + "/cylindrical_pincell.xml";

    shared_ptr<Constructive_Solid_Geometry> solid = get_solid_geometry(xml_input_filepath);

    shared_ptr<Distance> cartesian
        = make_shared<Cartesian_Distance>(dimension);
    shared_ptr<Distance> optical
        = make_shared<Optical_Distance>(dimension,
                                        solid);

    double const st_fuel = 1.0;
    double const st_mod = 2.0;

    {
        int test_case = 1;
        vector<double> const initial = {0, 0};
        vector<double> const final = {2, 1};

        double dist = sqrt(5);
        vector<double> grad = {2. / sqrt(5.), 1. / sqrt(5.)};
        vector<double> double_grad = {1. / (5 * sqrt(5.)), -2. / (5. * sqrt(5.)), -2. / (5. * sqrt(5.)), 4. / (5. * sqrt(5.))};
        
        double opt_dist = 1. * st_fuel + (sqrt(5.) - 1.) * st_mod;
        vector<double> opt_grad = vf::multiply(grad,
                                                st_mod);
        vector<double> opt_double_grad = vf::multiply(double_grad,
                                                       st_mod);
        
        checksum += test_distance(test_case,
                                  cartesian,
                                  dist,
                                  grad,
                                  double_grad,
                                  initial,
                                  final);
        
        checksum += test_distance(test_case,
                                  optical,
                                  opt_dist,
                                  opt_grad,
                                  opt_double_grad,
                                  initial,
                                  final);
    }
    
    return checksum;
}
