#include "Check_Equality.hh"
#include "Random_Number_Generator.hh"
#include "Trilinos_Dense_Solve.hh"

namespace ce = Check_Equality;

using namespace std;

int test_trilinos_dense(int number_of_points)
{
    int checksum = 0;
    
    Random_Number_Generator<double> rng(-1,
                                        1);
    Trilinos_Dense_Solve solver;

    vector<double> const a_data = rng.vector(number_of_points * number_of_points);
    vector<double> const b_data = rng.vector(number_of_points);

    vector<double> x_epetra(number_of_points);
    vector<double> x_amesos(number_of_points);

    {
        vector<double> a(a_data);
        vector<double> b(b_data);
        
        solver.epetra_solve(a,
                            b,
                            x_epetra,
                            number_of_points);
    }

    {
        vector<double> a(a_data);
        vector<double> b(b_data);
        
        solver.amesos_dense_solve(a,
                                  b,
                                  x_amesos,
                                  number_of_points);
    }

    if (!ce::approx(x_epetra, x_amesos, 1e-10))
    {
        checksum += 1;
    }

    return checksum;
}

int main()
{
    int checksum = 0;

    checksum += test_trilinos_dense(10);
    checksum += test_trilinos_dense(100);
    checksum += test_trilinos_dense(1000);
    
    return checksum;
}
