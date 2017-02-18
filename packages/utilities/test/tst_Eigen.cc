#include <iomanip>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "Linear_Algebra.hh"
#include "Random_Number_Generator.hh"
#include "Timer.hh"
#include "Trilinos_Dense_Solve.hh"

using namespace std;

#define global_size 2
#define num_tests 1000

#if global_size < 200
typedef Eigen::Map<Eigen::Matrix<double, global_size, global_size, Eigen::RowMajor> const> EMatrixRC;
typedef Eigen::Map<Eigen::Matrix<double, global_size, 1> const> EVectorC;
typedef Eigen::Map<Eigen::Matrix<double, global_size, 1> > EVector;
typedef Eigen::FullPivLU<Eigen::Matrix<double, global_size, global_size, Eigen::RowMajor> > ELU;

void solve(vector<double> const &a_data,
           vector<double> const &b_data,
           vector<double> &x_data)
{
    EMatrixRC A(&a_data[0]);
    EVectorC b(&b_data[0]);
    EVector x(&x_data[0]);
    x = A.fullPivLu().solve(b);
}

template<int const size>
class Solver
{
public:
    
    Solver()
    {
    }
    
    void solve(vector<double> const &a_data,
               vector<double> const &b_data,
               vector<double> &x_data)
    {
        typedef Eigen::Map<Eigen::Matrix<double, size, size, Eigen::RowMajor> const> EMatrixRC2;
        typedef Eigen::Map<Eigen::Matrix<double, size, 1> const> EVectorC2;
        typedef Eigen::Map<Eigen::Matrix<double, size, 1> > EVector2;
        
        EMatrixRC A(&a_data[0]);
        EVectorC b(&b_data[0]);
        EVector x(&x_data[0]);
        x = A.fullPivLu().solve(b);
    }
};
#endif

void solve2(vector<double> const &a_data,
           vector<double> const &b_data,
           vector<double> &x_data)
{
    int const size = b_data.size();
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const> A(&a_data[0], size, size);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1> const> b(&b_data[0], size);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1> > x(&x_data[0], size);
    x = A.fullPivLu().solve(b);
}

int main()
{
    Random_Number_Generator<double> rng(-1, 1);
    vector<double> a = rng.vector(global_size * global_size);
    vector<double> b = rng.vector(global_size);
    vector<double> x(global_size);

    Timer timer;
    
    int w1 = 20;
    int w2 = 16;

#if global_size < 200
    timer.start();
    for (int i = 0; i < num_tests; ++i)
    {
        solve(a, b, x);
    }
    timer.stop();
    cout << setw(w1) << "Eigen, fixed";
    cout << setw(w2) << timer.time();
    cout << endl;

    int const class_size = ;
    Solver<class_size> solver;
    timer.start();
    for (int i = 0; i < num_tests; ++i)
    {
        solver.solve(a, b, x);
    }
    timer.stop();
    cout << setw(w1) << "Eigen, class";
    cout << setw(w2) << timer.time();
    cout << endl;
    
    // {
    //     EMatrixRC A(&a[0]);
    //     EVectorC B(&b[0]);
    //     EVector X(&x[0]);
    //     ELU lu = A.fullPivLu();
    //     timer.start();
    //     for (int i = 0; i < num_tests; ++i)
    //     {
    //         X = lu.solve(B);
    //     }
    //     timer.stop();
    //     cout << setw(w1) << "Eigen, stored";
    //     cout << setw(w2) << timer.time();
    //     cout << endl;
    // }
#endif
    
    timer.start();
    for (int i = 0; i < num_tests; ++i)
    {
        solve2(a, b, x);
    }
    timer.stop();
    cout << setw(w1) << "Eigen, Dynamic";
    cout << setw(w2) << timer.time();
    cout << endl;

    if (global_size < 6)
    {
        timer.start();
        for (int i = 0; i < num_tests; ++i)
        {
            Linear_Algebra::linear_solve(a, b, x);
        }
        timer.stop();
        cout << setw(w1) << "Analytic";
        cout << setw(w2) << timer.time();
        cout << endl;
    }
    
    Trilinos_Dense_Solve tds;
    timer.start();
    for (int i = 0; i < num_tests; ++i)
    {
        tds.epetra_solve(a, b, x, global_size);
    }
    timer.stop();
    cout << setw(w1) << "Epetra";
    cout << setw(w2) << timer.time();
    cout << endl;
}
