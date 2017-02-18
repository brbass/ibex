#include "Dense_Solver_Factory.hh"

#include "Direct_Dense_Solver.hh"
#include "Eigen_Dense_Solver.hh"
#include "Eigen_Fixed_Dense_Solver.hh"
#include "Trilinos_Dense_Solver.hh"

// Explicitly create fixed dense solvers for later use
typedef Eigen_Fixed_Dense_Solver<1, Dense_Solver_Factory::Scalar> EFDS1;
typedef Eigen_Fixed_Dense_Solver<2, Dense_Solver_Factory::Scalar> EFDS2;
typedef Eigen_Fixed_Dense_Solver<3, Dense_Solver_Factory::Scalar> EFDS3;
typedef Eigen_Fixed_Dense_Solver<4, Dense_Solver_Factory::Scalar> EFDS4;
typedef Eigen_Fixed_Dense_Solver<5, Dense_Solver_Factory::Scalar> EFDS5;
typedef Eigen_Fixed_Dense_Solver<6, Dense_Solver_Factory::Scalar> EFDS6;
typedef Eigen_Fixed_Dense_Solver<7, Dense_Solver_Factory::Scalar> EFDS7;
typedef Eigen_Fixed_Dense_Solver<8, Dense_Solver_Factory::Scalar> EFDS8;
typedef Eigen_Fixed_Dense_Solver<9, Dense_Solver_Factory::Scalar> EFDS9;
typedef Eigen_Fixed_Dense_Solver<10, Dense_Solver_Factory::Scalar> EFDS10;
typedef Eigen_Fixed_Dense_Solver<11, Dense_Solver_Factory::Scalar> EFDS11;
typedef Eigen_Fixed_Dense_Solver<12, Dense_Solver_Factory::Scalar> EFDS12;
typedef Eigen_Fixed_Dense_Solver<13, Dense_Solver_Factory::Scalar> EFDS13;
typedef Eigen_Fixed_Dense_Solver<14, Dense_Solver_Factory::Scalar> EFDS14;
typedef Eigen_Fixed_Dense_Solver<15, Dense_Solver_Factory::Scalar> EFDS15;
typedef Eigen_Fixed_Dense_Solver<16, Dense_Solver_Factory::Scalar> EFDS16;
typedef Eigen_Fixed_Dense_Solver<17, Dense_Solver_Factory::Scalar> EFDS17;
typedef Eigen_Fixed_Dense_Solver<18, Dense_Solver_Factory::Scalar> EFDS18;
typedef Eigen_Fixed_Dense_Solver<19, Dense_Solver_Factory::Scalar> EFDS19;
typedef Eigen_Fixed_Dense_Solver<20, Dense_Solver_Factory::Scalar> EFDS20;

using namespace std;

Dense_Solver_Factory::
Dense_Solver_Factory()
{
}

shared_ptr<Dense_Solver<Dense_Solver_Factory::Scalar> > Dense_Solver_Factory::
get_solver(int size,
           Type type) const
{
    switch (type)
    {
    case Type::DEFAULT:
        return get_default(size);
    case Type::DIRECT:
        return make_shared<Direct_Dense_Solver<double> >(size);
    case Type::EIGEN:
        return make_shared<Eigen_Dense_Solver<double> >(size);
    case Type::EIGEN_FIXED:
        return get_fixed_eigen_solver(size);
    case Type::EPETRA:
        return make_shared<Trilinos_Dense_Solver>(size);
    }
}

shared_ptr<Dense_Solver<Dense_Solver_Factory::Scalar> > Dense_Solver_Factory::
get_default(int size) const
{
    if (size <= 5)
    {
        return get_solver(size,
                          Type::DIRECT);
    }
    else if (size <= 20)
    {
        return get_solver(size,
                          Type::EIGEN_FIXED);
    }
    else
    {
        return get_solver(size,
                          Type::EPETRA);
    }
}

shared_ptr<Dense_Solver<Dense_Solver_Factory::Scalar> > Dense_Solver_Factory::
get_fixed_eigen_solver(int size) const
{
    switch (size)
    {
    case 1:
        return make_shared<EFDS1>();
    case 2:
        return make_shared<EFDS2>();
    case 3:
        return make_shared<EFDS3>();
    case 4:
        return make_shared<EFDS4>();
    case 5:
        return make_shared<EFDS5>();
    case 6:
        return make_shared<EFDS6>();
    case 7:
        return make_shared<EFDS7>();
    case 8:
        return make_shared<EFDS8>();
    case 9:
        return make_shared<EFDS9>();
    case 10:
        return make_shared<EFDS10>();
    case 11:
        return make_shared<EFDS11>();
    case 12:
        return make_shared<EFDS12>();
    case 13:
        return make_shared<EFDS13>();
    case 14:
        return make_shared<EFDS14>();
    case 15:
        return make_shared<EFDS15>();
    case 16:
        return make_shared<EFDS16>();
    case 17:
        return make_shared<EFDS17>();
    case 18:
        return make_shared<EFDS18>();
    case 19:
        return make_shared<EFDS19>();
    case 20:
        return make_shared<EFDS20>();            
    }
}
