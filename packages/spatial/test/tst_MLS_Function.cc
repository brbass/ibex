#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

#include "Cartesian_Distance.hh"
#include "Check.hh"
#include "Check_Equality.hh"
#include "Compact_Gaussian_RBF.hh"
#include "Linear_MLS_Function.hh"
#include "Meshless_Function_Factory.hh"
#include "Random_Number_Generator.hh"
#include "RBF_Function.hh"
#include "Timer.hh"

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
                 return (std::abs(func1->position()[0] - position)
                         < std::abs(func2->position()[0] - position));
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
        double tolerance = 1e-9;
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
        double tolerance = 1e-8;
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
            double error = std::abs(test_value);
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

// int test_values()
// {
//     int checksum = 0;

//     Meshless_Function_Factory meshless_factory;
//     int dimension = 1;
//     int order = 1;
//     double radius_num_intervals = 5.1;
//     vector<int> dimensional_points(dimension, 11);
//     vector<vector<double> > limits(dimension);
//     limits[0] = {0, 2};
//     string rbf_type = "wendland11";
//     vector<shared_ptr<Meshless_Function> > functions;
//     meshless_factory.get_cartesian_mls_functions(order,
//                                                  dimension,
//                                                  radius_num_intervals,
//                                                  dimensional_points,
//                                                  limits,
//                                                  rbf_type,
//                                                  functions);
//     shared_ptr<Meshless_Function> function = functions[0];
//     cout << function->value({0.1}) << endl;
//     cout << function->d_value(0, {0.1}) << endl;
    
//     return checksum;
// }

int test_all_values(int order,
                    int dimension)
{
    int checksum = 0;

    cout << "dimension: " << dimension << "\t" << "order: " << order << endl;
    
    // Set preliminary data
    double tolerance = order == 1 ? 1e-9 : 1e-8;
    double grad_tolerance = order == 1 ? 1e-9 : 2e-8;
    double radius_num_intervals = 3.5;
    string rbf_type = "wendland11";
    
    // Set limits and number of points
    vector<vector<double> > limits(dimension, vector<double>(2));
    for (int d = 0; d < dimension; ++d)
    {
        limits[d][0] = -1;
        limits[d][1] = 1;
    }
    vector<int> dimensional_points(dimension, 20);
    
    // Get MLS functions
    vector<shared_ptr<Meshless_Function> > functions;
    Meshless_Function_Factory factory;
    factory.get_cartesian_mls_functions(order,
                                        dimension,
                                        radius_num_intervals,
                                        dimensional_points,
                                        limits,
                                        rbf_type,
                                        functions);
    
    // Test values of two methods
    Random_Number_Generator<double> position_rng(-1, 1, 295 /*seed*/);
    Random_Number_Generator<int> function_rng(0, functions.size() - 1, 104 /*seed*/);
    int num_tests = 1000;
    int num_inside = 0;
    Timer timer;
    double new_time = 0;
    double old_time = 0;
    for (int t = 0; t < num_tests; ++t)
    {
        // Get random point and function
        shared_ptr<Meshless_Function> func = functions[function_rng.scalar()];
        vector<double> position = position_rng.vector(dimension);
        
        if (func->inside_radius(position))
        {
            num_inside += 1;

            // Get vector of values
            vector<int> indices;
            vector<double> vals;
            vector<vector<double> > grad_vals;
            timer.start();
            func->gradient_values(position,
                                  indices,
                                  vals,
                                  grad_vals);
            timer.stop();
            new_time += timer.time();
            
            // Test against individual values
            for (int i = 0; i < indices.size(); ++i)
            {
                timer.start();
                int j = indices[i];
                double val = functions[j]->value(position);
                vector<double> grad_val = functions[j]->gradient_value(position);
                timer.stop();
                old_time += timer.time();
                
                if (!ce::approx(val, vals[i], tolerance))
                {
                    cout << "MLS value failed for test " << t << " and value " << i << endl;
                    cout << "dimension: " << dimension << "\torder: " << order << endl;
                    cout << "\terror: " << val - vals[i] << endl;
                    checksum += 1;
                }
                if (!ce::approx(grad_val, grad_vals[i], grad_tolerance))
                {
                    cout << "MLS grad value failed for test " << t << " and value " << i << endl;
                    cout << "\terror: ";
                    for (int d = 0; d < dimension; ++d)
                    {
                        cout << grad_val[d] - grad_vals[i][d] << "\t";
                    }
                    cout << endl;
                    checksum += 1;
                }
            }
        }
    }
    Assert(num_inside > 10);

    cout << "\told method time: " << old_time << endl;
    cout << "\tnew method time: " << new_time << endl;
    
    return checksum;
}

int main()
{
    int checksum = 0;

    // checksum += test_1d();
    for (int order = 1; order <=2; ++ order)
    {
        for (int d = 1; d <= 3; ++d)
        {
            checksum += test_all_values(order,
                                        d);
        }
    }
    
    return checksum;
}
