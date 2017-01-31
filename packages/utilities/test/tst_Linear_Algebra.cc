#include <iostream>

#include "Check_Equality.hh"
#include "Linear_Algebra.hh"
#include "Random_Number_Generator.hh"
#include "Timer.hh"
#include "Trilinos_Dense_Solve.hh"

namespace ce = Check_Equality;
namespace la = Linear_Algebra;

using namespace std;

namespace // anonymous
{
    Random_Number_Generator<double> rng(0,    // min
                                        1,    // max
                                        948); // seed
}

int test_against_trilinos(int size)
{
    int checksum = 0;

    Trilinos_Dense_Solve trilinos_solver;

    int num_tests = 10;
    double tolerance = 1e-12;

    for (int i = 0; i < num_tests; ++i)
    {
        // Get random matrix and LHS
        vector<double> const a = rng.vector(size * size);
        vector<double> const b = rng.vector(size);
        
        // Solve using Linear_Algebra
        vector<double> x_ls(size);
        la::linear_solve(a, b, x_ls);

        // Solve using Trilinos
        vector<double> x_tr(size);
        vector<double> a_tr = a;
        vector<double> b_tr = b;
        trilinos_solver.epetra_solve(a_tr, b_tr, x_tr, size);

        // Compare results
        if (!(ce::approx(x_tr, x_ls, tolerance)))
        {
            cerr << "linear solve incorrect for size " << size << endl;
            cerr << "\texpected: " << x_tr[0] << "\tcalculated: " << x_ls[0] << endl;
            checksum += 1;
        }
    }

    return checksum;
}

int main()
{
    int checksum = 0;

    // Check linear solve for matrices of size 1 to 6
    for (int i = 1; i < 6; ++i)
    {
        checksum += test_against_trilinos(i);
    }
    
    return checksum;
}
