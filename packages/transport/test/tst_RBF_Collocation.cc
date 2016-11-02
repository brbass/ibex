#include <iomanip>
#include <iostream>
#include <memory>

#include <mpi.h>

#include "Angular_Discretization_Parser.hh"
#include "Boundary_Source_Parser.hh"
#include "Check.hh"
#include "Check_Equality.hh"
#include "Energy_Discretization_Parser.hh"
#include "Material_Parser.hh"
#include "Random_Number_Generator.hh"
#include "RBF_Collocation_Sweep.hh"
#include "Solver_Parser.hh"
#include "Spatial_Discretization_Parser.hh"
#include "Sweep_Parser.hh"
#include "Transport_Discretization.hh"
#include "Vector_Operator.hh"

using namespace std;
namespace ce = Check_Equality;

shared_ptr<RBF_Collocation_Sweep> get_rbf_sweep(string xml_input_filename)
{
    pugi::xml_document input_document;
    
    if (!input_document.load_file(xml_input_filename.c_str()))
    {
        cout << "Could not open xml input file \"" + xml_input_filename + "\"" << endl;
        return shared_ptr<RBF_Collocation_Sweep>();
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
    shared_ptr<RBF_Discretization> rbf
        = dynamic_pointer_cast<RBF_Discretization>(spatial);
    Assert(rbf);
    
    shared_ptr<Transport_Discretization> transport
        = make_shared<Transport_Discretization>(spatial,
                                                angular,
                                                energy);
    
    Sweep_Parser sweep_parser(input_file,
                              spatial,
                              angular,
                              energy,
                              transport);
    
    shared_ptr<Sweep_Operator> sweep_operator = sweep_parser.get_ptr();
    shared_ptr<RBF_Collocation_Sweep> rbf_collocation
        = dynamic_pointer_cast<RBF_Collocation_Sweep>(sweep_operator);
    Assert(rbf_collocation);

    return rbf_collocation;
}

int test_free_streaming(string input_folder)
{
    int checksum = 0;
    
    string input_file = input_folder + "/square_free_streaming.xml";
    
    shared_ptr<RBF_Collocation_Sweep> sweep = get_rbf_sweep(input_file);
    shared_ptr<RBF_Discretization> spatial = sweep->rbf_discretization();
    shared_ptr<Angular_Discretization> angular = sweep->angular_discretization();
    shared_ptr<Energy_Discretization> energy = sweep->energy_discretization();
    shared_ptr<Transport_Discretization> transport = sweep->transport_discretization();
    
    sweep->set_include_boundary_source(true);
    
    int psi_size = transport->psi_size();
    int number_of_augments = transport->number_of_augments();
    
    vector<double> x(psi_size + number_of_augments, 0);
    
    (*sweep)(x);
    
    x.resize(psi_size);
    vector<double> solution(psi_size, 1.0);
    
    if (!ce::approx(x, solution, 1e-6))
    {
        checksum += 1;
        cout << "free streaming failed" << endl;

        for (int i = 0; i < 100; ++i)
        {
            cout << x[i] << "\t";
        }
        cout << endl;
    }
    
    return checksum;
}

int test_manufactured(string input_file)
{
    int checksum = 0;
    
    shared_ptr<RBF_Collocation_Sweep> sweep = get_rbf_sweep(input_file);
    shared_ptr<RBF_Discretization> spatial = sweep->rbf_discretization();
    shared_ptr<Angular_Discretization> angular = sweep->angular_discretization();
    shared_ptr<Energy_Discretization> energy = sweep->energy_discretization();
    shared_ptr<Transport_Discretization> transport = sweep->transport_discretization();
    
    sweep->set_include_boundary_source(false);
    
    int psi_size = transport->psi_size();
    int number_of_augments = transport->number_of_augments();
    int dimension = spatial->dimension();
    int number_of_points = spatial->number_of_points();
    int number_of_groups = energy->number_of_groups();
    int number_of_ordinates = angular->number_of_ordinates();
    vector<double> wavenumber(dimension, 2 * M_PI);
    
    vector<double> source(psi_size, 0);
    vector<double> solution(psi_size, 0.0);
    
    for (int i = 0; i < number_of_points; ++i)
    {
        shared_ptr<RBF_Point> point = spatial->rbf_point(i);
        vector<double> const position = point->position();
        vector<double> const sigma_t = point->material()->sigma_t();
        
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            vector<double> const direction = angular->direction(o);
            
            double streaming = 0;
            for (int d1 = 0; d1 < dimension; ++d1)
            {
                double product = wavenumber[d1] * direction[d1] * cos(wavenumber[d1] * position[d1]);
                
                for (int d2 = 0; d2 < dimension; ++d2)
                {
                    if (d1 != d2)
                    {
                        product *= sin(wavenumber[d2] * position[d2]);
                    }
                }
                
                streaming += product;
            }
            
            double absorption = 1;
            for (int d = 0; d < dimension; ++d)
            {
                absorption *= sin(wavenumber[d] * position[d]);
            }
            
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k = g + number_of_groups * (o + number_of_ordinates * i);
                
                source[k] = streaming + sigma_t[g] * absorption;
                
                solution[k] = absorption;
            }
        }
    }
    
    vector<double> x(source);
    x.resize(psi_size + number_of_augments, 0);
    (*sweep)(x);
    x.resize(psi_size);

    double norm = 0;
    if (!ce::norm_approx(x, solution, 1e-4, norm))
    {
        checksum += 1;

        Random_Number_Generator<int> point(0, number_of_points, 0);
        Random_Number_Generator<int> ordinate(0, number_of_ordinates, 1);
        Random_Number_Generator<int> group(0, number_of_groups, 2);
        
        cout << "manufactured solution failed for ";
        cout << input_file;
        cout << endl;
        cout << "\tnorm:  " << norm << endl;
        
        int num_to_print = 100;
        int ws = 6;
        int w = 16;
        cout << setw(ws) << "i";
        cout << setw(ws) << "g";
        cout << setw(ws) << "o";
        cout << setw(w) << "source";
        cout << setw(w) << "solution";
        cout << setw(w) << "expected";
        cout << endl;
        
        for (int n = 0; n < num_to_print; ++n)
        {
            int i = point.scalar();
            int o = ordinate.scalar();
            int g = group.scalar();
            int k = g + number_of_groups * (o + number_of_ordinates * i);

            cout << setw(ws) << i;
            cout << setw(ws) << g;
            cout << setw(ws) << o;
            cout << setw(w) << source[k];
            cout << setw(w) << x[k];
            cout << setw(w) << solution[k];
            cout << endl;
        }
    }
    
    return checksum;
}

int main(int argc, char **argv)
{
    int checksum = 0;

    MPI_Init(&argc, &argv);

    if (argc != 2)
    {
        cerr << "usage: tst_RBF_Collocation [input_folder]" << endl;
        return 1;
    }

    string input_folder = argv[1];
    
    // checksum += test_free_streaming(input_folder);
    // checksum += test_pincell_manufactured(input_folder, "cartesian");
    // checksum += test_manufactured(input_folder + "/pincell_cartesian.xml");
    // checksum += test_manufactured(input_folder + "/pincell_optical.xml");
    // checksum += test_manufactured(input_folder + "/split_cartesian.xml");
    checksum += test_manufactured(input_folder + "/split_optical.xml");
    
    MPI_Finalize();
    
    return checksum;
}
