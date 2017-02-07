#include <algorithm>
#include <memory>
#include <vector>

#include "Cartesian_Distance.hh"
#include "Compact_Gaussian_RBF.hh"
#include "Linear_MLS_Function.hh"
#include "RBF_Function.hh"

int test_1d()
{
    int dimension = 1;
    int number_of_points = 21;
    double dx = 0.2;
    int number_of_neighbors = 4;
    double radius = dx * number_of_neighbors;
    
    shared_ptr<RBF> rbf
        = make_shared<Compact_Gaussian_RBF>();
    shared_ptr<Distance> distance
        = make_shared<Cartesian_Distance(dimension);
    
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
    
    vector<shared_ptr<Meshless_Function> > mls_functions(number_of_points);
    for (int i = 0; i < number_of_points; ++i)
    {
        int starting_j = i - number_of_neighbors;
        int ending_j = i + number_of_neighbors;

        double position = dx * i;
        
        starting_j = starting_j < 0 ? 0 : starting_j;
        ending_j = ending_j > number_of_points - 1 ? number_of_points - 1 : ending_j;
        
        vector<shared_ptr<Meshless_Function> > local_rbf_functions;        
        for (int j = starting_j; j <= ending_j; ++j)
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
        
        for (shared_ptr<Meshless_Function> local_function : local_rbf_functions)
        {
            cout << local_rbf_functions->position()[0] << endl;
        }
        
        
    }
}

int main()
{
    int checksum = 0;

    checksum += test_1d();

    return checksum;
}
