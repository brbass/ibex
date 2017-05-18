#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

#include "Cartesian_Distance.hh"
#include "Check_Equality.hh"
#include "Compact_Gaussian_RBF.hh"
#include "Linear_MLS_Function.hh"
#include "Random_Number_Generator.hh"
#include "RBF_Function.hh"
#include "Weak_Spatial_Discretization_Factory.hh"

using namespace std;
namespace ce = Check_Equality;

int test_1d()
{
    int checksum = 0;

    // Set up physical data
    int dimension = 1;
    int number_of_points = 21;
    double dx = 0.2;
    int number_of_neighbors = 6;
    double xmin = 0.;
    double xmax = dx * (number_of_points - 1);

    // Radius is halved so that a basis function 2r away will be the final neighbor
    double radius = 0.5 * dx * number_of_neighbors;
    
    // Create RBF and distance
    shared_ptr<RBF> rbf
        = make_shared<Compact_Gaussian_RBF>();
    shared_ptr<Distance> distance
        = make_shared<Cartesian_Distance>(dimension);

    // Create weighting functions
    vector<shared_ptr<Meshless_Function> > rbf_functions(number_of_points);
    for (int i = 0; i < number_of_points; ++i)
    {
        double shape = rbf->radius() / radius;
        vector<double> position(1, dx * i);
        
        rbf_functions[i]
            = make_shared<RBF_Function>(i,
                                        shape,
                                        position,
                                        rbf,
                                        distance);
    }

    // Create MLS functions
    vector<shared_ptr<Meshless_Function> > mls_functions(number_of_points);
    for (int i = 0; i < number_of_points; ++i)
    {
        double position = dx * i;

        // Get starting and ending neighbor indices
        int starting_j = i - number_of_neighbors;
        int ending_j = i + number_of_neighbors;

        starting_j = starting_j < 0 ? 0 : starting_j;
        ending_j = ending_j > number_of_points - 1 ? number_of_points - 1 : ending_j;
        
        // Get local functions that overlap and sort by distance to center
        vector<shared_ptr<Meshless_Function> > local_rbf_functions;             for (int j = starting_j; j <= ending_j; ++j)
        {
            local_rbf_functions.push_back(rbf_functions[j]);
        }
        
        sort(local_rbf_functions.begin(),
             local_rbf_functions.end(),
             [position](shared_ptr<Meshless_Function> func1,
                        shared_ptr<Meshless_Function> func2)
             {
                 return (abs(func1->position()[0] - position)
                         < abs(func2->position()[0] - position));
             });

        // Create MLS function
        mls_functions[i]
            = make_shared<Linear_MLS_Function>(local_rbf_functions);
    }

    // General test data
    int num_tests = 1000;
    Random_Number_Generator<double> rng(xmin, xmax, 294 /*seed*/);

    // Test that sum of values is one everywhere
    {
        bool failed = false;
        double tolerance = 1.e-10;
        double largest_error = 0.;
        for (int t = 0; t < num_tests; ++t)
        {
            vector<double> test_position(dimension, rng.scalar());
        
            double test_value = 0.;
            for (int i = 0; i < number_of_points; ++i)
            {
                double value = mls_functions[i]->value(test_position);
                test_value += value;
            }
            double error = test_value - 1.;
            if (error > largest_error)
            {
                largest_error = error;
            }
            if (!ce::approx(error, 0., tolerance))
            {
                checksum += 1;
                if (!failed)
                {
                    failed = true;
                    cout << "MLS value error for test " << t << ": " << error << endl;
                }
            }
        }
        cout << "largest error for MLS value: " << largest_error << endl;
    }
    
    // Test that sum of derivatives is zero everywhere
    {
        bool failed = false;
        double tolerance = 1.e-8;
        double largest_error = 0.;
        for (int t = 0; t < num_tests; ++t)
        {
            vector<double> test_position(dimension, rng.scalar());
        
            double test_value = 0.;
            for (int i = 0; i < number_of_points; ++i)
            {
                double value = mls_functions[i]->d_value(0,
                                                         test_position);
                test_value += value;
            }
            double error = abs(test_value);
            if (error > largest_error)
            {
                largest_error = error;
            }
            if (!ce::approx(test_value, 0., tolerance))
            {
                checksum += 1;
                if (!failed)
                {
                    failed = true;
                    cout << "MLS d_value error for test " << t << ": " << test_value << endl;
                }
            }
        }
        cout << "largest error for MLS d_value: " << largest_error << endl;
    } 
   
    return checksum;
}

int test_values()
{
    int checksum = 0;

    Meshless_Function_Factory meshless_factory;
    int dimension = 1;
    double radius_num_intervals = 5.1;
    vector<int> dimensional_points(dimension, 11);
    vector<vector<double> > limits(dimension);
    limits[0] = {0, 2};
    string rbf_type = "wendland11";
    vector<shared_ptr<Meshless_Function> > functions;
    meshless_factory.get_cartesian_mls_functions(dimension,
                                                 radius_num_intervals,
                                                 dimensional_points,
                                                 limits,
                                                 rbf_type,
                                                 functions);
    shared_ptr<Meshless_Function> function = functions[0];
    cout << function->value({0.1}) << endl;
    cout << function->d_value(0, {0.1}) << endl;
    
    return checksum;
}

int main()
{
    int checksum = 0;

    checksum += test_1d();
    
    return checksum;
}
