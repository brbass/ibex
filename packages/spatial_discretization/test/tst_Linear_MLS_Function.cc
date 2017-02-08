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

using namespace std;
namespace ce = Check_Equality;

int test_1d()
{
    int checksum = 0;

    // Set up physical data
    int dimension = 1;
    int number_of_points = 21;
    double dx = 0.2;
    int number_of_neighbors = 4;
    double radius = dx * number_of_neighbors;
    double xmin = 0.;
    double xmax = dx * (number_of_points - 1);

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
            = make_shared<RBF_Function>(shape,
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
        double tolerance = 1.e-6;
        for (int t = 0; t < num_tests; ++t)
        {
            vector<double> test_position(dimension, rng.scalar());
        
            double test_value = 0.;
            for (int i = 0; i < number_of_points; ++i)
            {
                double value = mls_functions[i]->value(test_position);
                test_value += value;
            }
            if (!ce::approx(test_value, 1., tolerance))
            {
                checksum += 1;
                if (!failed)
                {
                    failed = true;
                    cout << "MLS value error: " << test_value - 1. << endl;
                }
            }
        }
    }
    
    // Test that sum of derivatives is zero everywhere
    {
        bool failed = false;
        double tolerance = 2.e-5;
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
            if (!ce::approx(test_value, 0., tolerance))
            {
                checksum += 1;
                if (!failed)
                {
                    failed = true;
                    cout << "MLS d_value error: " << test_value << endl;
                }
            }
        }
    } 
   
    return checksum;
}

int main()
{
    int checksum = 0;

    checksum += test_1d();
    
    return checksum;
}
