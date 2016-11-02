#ifndef Trilinos_Dense_Solve_hh
#define Trilinos_Dense_Solve_hh

#include <vector>

/*
  Dense matrix solution routines from GNU Scientific Library
  
  Not optimized for matrix reuse - does not store factorizations
*/
class Trilinos_Dense_Solve
{
private:
    
public:

    // Constructor
    Trilinos_Dense_Solve();

    // Solve using Epetra_SerialDenseSolver (LU decomposition)
    void epetra_solve(std::vector<double> &a_data,
                      std::vector<double> &b_data,
                      std::vector<double> &x_data,
                      unsigned number_of_elements);

    // Solve using Amesos Klu
    // Not optimal, as it uses a sparse matrix representation for a dense problem
    void amesos_dense_solve(std::vector<double> &a_data,
                            std::vector<double> &b_data,
                            std::vector<double> &x_data,
                            unsigned number_of_elements);
    
    // Solve using Aztec GMRES
    // Not optimal, as it uses a sparse matrix representation for a dense problem
    void aztec_dense_solve(std::vector<double> &a_data,
                           std::vector<double> &b_data,
                           std::vector<double> &x_data,
                           unsigned number_of_elements);
};

#endif
