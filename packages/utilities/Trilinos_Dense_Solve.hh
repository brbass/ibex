#ifndef Trilinos_Dense_Solve_hh
#define Trilinos_Dense_Solve_hh

#include "Dense_Solve.hh"

/*
  Three methods from Trilinos for matrix solution
*/
class Trilinos_Dense_Solve : public Dense_Solve<double>
{
public:

    enum class Solver
    {
        EPETRA,
            AMESOS,
            AZTEC
            };
    
    Trilinos_Dense_Solve(int size,
                         Solver solver = Solver::EPETRA);
    
    virtual void solve(std::vector<double> &a_data,
                       std::vector<double> &b_data,
                       std::vector<double> &x_data) const override;

    virtual int size() const override
    {
        return size_;
    }

    // Solve using Epetra_SerialDenseSolver (LU decomposition)
    void epetra_solve(std::vector<double> &a_data,
                      std::vector<double> &b_data,
                      std::vector<double> &x_data) const;
    
    // Solve using Amesos Klu
    // Not optimal, as it uses a sparse matrix representation for a dense problem
    void amesos_solve(std::vector<double> &a_data,
                      std::vector<double> &b_data,
                      std::vector<double> &x_data) const;
    
    // Solve using Aztec GMRES
    // Not optimal, as it uses a sparse matrix representation for a dense problem
    void aztec_solve(std::vector<double> &a_data,
                     std::vector<double> &b_data,
                     std::vector<double> &x_data) const;

private:

    int size_;
    Solver solver_;
};

#endif
