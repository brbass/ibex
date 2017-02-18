#include <mpi.h>

#include "Check_Equality.hh"
#include "Random_Number_Generator.hh"
#include "Trilinos_Dense_Solve.hh"

namespace ce = Check_Equality;

using namespace std;

int test_trilinos_dense(int number_of_points)
{
    int checksum = 0;
    
    Random_Number_Generator<double> rng(-1,
                                        1,
                                        402 * number_of_points);
    Trilinos_Dense_Solve solver(number_of_points);

    vector<double> const a_data = rng.vector(number_of_points * number_of_points);
    vector<double> const b_data = rng.vector(number_of_points);

    vector<double> x_epetra(number_of_points);
    vector<double> x_amesos(number_of_points);

    {
        vector<double> a(a_data);
        vector<double> b(b_data);
        
        solver.epetra_solve(a,
                            b,
                            x_epetra);
    }

    {
        vector<double> a(a_data);
        vector<double> b(b_data);
        
        solver.amesos_solve(a,
                            b,
                            x_amesos);
    }

    if (!ce::approx(x_epetra, x_amesos, 1e-9))
    {
        cerr << "tst_Trilinos: results differ from amesos to epetra";
        cerr << endl;
        cerr << "error of first two points: ";
        cerr << x_epetra[0] - x_amesos[0];
        cerr << ", ";
        cerr << x_epetra[1] - x_amesos[1];
        cerr << endl;
            
        checksum += 1;
    }

    return checksum;
}

int main(int argc, char **argv)
{
    int checksum = 0;

    MPI_Init(&argc, &argv);
    
    checksum += test_trilinos_dense(10);
    checksum += test_trilinos_dense(100);
    checksum += test_trilinos_dense(1000);

    MPI_Finalize();
    
    return checksum;
}
