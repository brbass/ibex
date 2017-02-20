#include <algorithm>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <memory>
#include <string>
#include <vector>

#include "Check_Equality.hh"
#include "Dense_Solver_Factory.hh"
#include "Matrix_Functions.hh"
#include "Random_Number_Generator.hh"
#include "Timer.hh"

using namespace std;

Random_Number_Generator<double> rng(0, 1, 59123);
Timer timer;

struct Solver_Data
{
    Solver_Data(string description,
                shared_ptr<Dense_Solver<double> > solver):
        description(description),
        solver(solver)
    {
        time = 0;
        error = 0;
    }
    string description;
    double time;
    double error;
    double det;
    shared_ptr<Dense_Solver<double> > solver;
};

int test_linear_solve(vector<double> const &a_data,
                      vector<double> const &b_data,
                      vector<double> const &x_data,
                      Solver_Data &solver_data)
{
    int checksum = 0;

    double tolerance = 1e-8;
    
    vector<double> a_temp = a_data;
    vector<double> b_temp = b_data;
    vector<double> x_temp(x_data.size());
    
    timer.start();
    solver_data.solver->solve(a_temp,
                              b_temp,
                              x_temp);
    timer.stop();
    solver_data.time += timer.time();
    
    if (!Check_Equality::approx(x_temp, x_data, tolerance))
    {
        checksum += 1;
        cout << "solver " + solver_data.description + " failed to converge for size (" << x_data.size() << "); first element data:" << endl;
        cout << "\texpected: " << x_data[0] << "\tcalculated: " << x_temp[0] << "\terror:  " << x_temp[0] - x_data[0] << endl;
    }
    
    return checksum;
}

void get_linear_problem(int size,
                        vector<double> &a,
                        vector<double> &b,
                        vector<double> &x)
{
    a = rng.vector(size * size);
    x = rng.vector(size);
    b = Matrix_Functions::square_matrix_vector_product(size,
                                                       a,
                                                       x);
}

void get_solvers(int size,
                 vector<Solver_Data> &solvers,
                 bool include_trilinos = true)
{
    solvers.clear();
    
    Dense_Solver_Factory factory;
    
    if (size <= 5)
    {
        solvers.emplace_back("direct",
                             factory.get_solver(size,
                                                Dense_Solver_Factory::Type::DIRECT));
    }
    if (size <= 10)
    {
        solvers.emplace_back("eigen_fixed",
                             factory.get_solver(size,
                                                Dense_Solver_Factory::Type::EIGEN_FIXED));
    }
    solvers.emplace_back("eigen",
                         factory.get_solver(size,
                                            Dense_Solver_Factory::Type::EIGEN));
    if (include_trilinos)
    {
        solvers.emplace_back("epetra",
                             factory.get_solver(size,
                                                Dense_Solver_Factory::Type::EPETRA));
    }
}

int run_linear_solve(int number_of_tests,
                     vector<int> sizes)
{
    int checksum = 0;
    
    for (int size : sizes)
    {
        vector<Solver_Data> solvers;
        get_solvers(size,
                    solvers);
        
        for (int t = 0; t < number_of_tests; ++t)
        {
            vector<double> a;
            vector<double> b;
            vector<double> x;

            get_linear_problem(size,
                               a,
                               b,
                               x);
            
            for (Solver_Data &solver_data : solvers)
            {
                checksum += test_linear_solve(a, b, x, solver_data);
            }
        }

        sort(solvers.begin(), solvers.end(),
             [](const Solver_Data &a, const Solver_Data &b)
             {
                 return a.time < b.time;
             });

        cout << "solver times for n = ";
        cout << size;
        cout << endl;
        for (Solver_Data &solver : solvers)
        {
            int w = 16;
            cout << setw(w) << solver.description;
            cout << setw(w) << solver.time;
            cout << endl;
        }
    }
    
    return checksum;
}

int run_determinant()
{
    vector<Solver_Data> solvers;

    for (int size = 1; size <= 5; ++size)
    {
        vector<double> a = rng.vector(size * size);
        get_solvers(size,
                    solvers,
                    false);

        cout << "solver dets for n = ";
        cout << size;
        cout << endl;
        for (Solver_Data &solver : solvers)
        {
            solver.solver->initialize(a);
            solver.det = solver.solver->determinant();

            int w = 16;
            cout << setw(w) << solver.description;
            cout << setw(w) << solver.det;
            cout << endl;
        }

        for (Solver_Data &solver : solvers)
        {
            if (!Check_Equality::approx(solver.det, solvers[0].det, 1e-10))
            {
                cout << "determinants disagree for n = ";
                cout << size;
                cout << "\tdiff from first: ";
                cout << solver.det - solvers[0].det;
                cout << endl;
            }
        }
        
    }
}

int main(int argc, char **argv)
{
    int checksum = 0;
    
    MPI_Init(&argc, &argv);

    // Set to 1 for a quick test and 10-100 for more accurate timing
    int iteration_multiplier = 1;
    
    vector<int> small_sizes;
    for (int i = 1; i <= 20; ++i)
    {
        small_sizes.push_back(i);
    }
    
    checksum += run_linear_solve(100 * iteration_multiplier,
                                 small_sizes);
    
    vector<int> medium_sizes;
    for (int i = 30; i <= 200; i += 10)
    {
        medium_sizes.push_back(i);
    }
    checksum += run_linear_solve(10 * iteration_multiplier,
                                 medium_sizes);

    vector<int> large_sizes;
    for (int i = 300; i <= 1000; i += 100)
    {
        large_sizes.push_back(i);
    }
    checksum += run_linear_solve(1 * iteration_multiplier,
                                 large_sizes);
    
    checksum += run_determinant();
    
    MPI_Finalize();
    
    return checksum;
}
