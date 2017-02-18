#include <algorithm>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <memory>
#include <string>
#include <vector>

#include "Check_Equality.hh"
#include "Matrix_Functions.hh"
#include "Random_Number_Generator.hh"
#include "Timer.hh"

#include "Eigen_Fixed_LU_Decomposition.hh"

using namespace std;

Random_Number_Generator<double> rng(0, 1, 4120);
Timer timer;

typedef Eigen_Fixed_LU_Decomposition<2, double> LU_Test;


int main(int argc, char **argv)
{
    int checksum = 0;
    
    MPI_Init(&argc, &argv);

    vector<double> test_a = rng.vector(4);
    vector<double> test_b = rng.vector(2);
    vector<double> test_x(2);
    
    LU_Test lu_test(test_a);
    lu_test.solve(test_b,
                  test_x);

    cout << test_x[0] << endl;

    test_a[0] = 5;

    lu_test.solve(test_b,
                  test_x);
    
    MPI_Finalize();
    
    return checksum;
}
