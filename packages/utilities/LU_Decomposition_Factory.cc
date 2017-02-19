#include "LU_Decomposition_Factory.hh"

#include "Eigen_Fixed_LU_Decomposition.hh"
#include "Eigen_Fixed_Dense_LU_Decomposition.hh"
#include "Epetra_LU_Decomposition.hh"

// Explicitly create fixed dense solvers for later use
typedef Eigen_Fixed_LU_Decomposition<1, LU_Decomposition_Factory::Scalar> EFLU1;
typedef Eigen_Fixed_LU_Decomposition<2, LU_Decomposition_Factory::Scalar> EFLU2;
typedef Eigen_Fixed_LU_Decomposition<3, LU_Decomposition_Factory::Scalar> EFLU3;
typedef Eigen_Fixed_LU_Decomposition<4, LU_Decomposition_Factory::Scalar> EFLU4;
typedef Eigen_Fixed_LU_Decomposition<5, LU_Decomposition_Factory::Scalar> EFLU5;
typedef Eigen_Fixed_LU_Decomposition<6, LU_Decomposition_Factory::Scalar> EFLU6;
typedef Eigen_Fixed_LU_Decomposition<7, LU_Decomposition_Factory::Scalar> EFLU7;
typedef Eigen_Fixed_LU_Decomposition<8, LU_Decomposition_Factory::Scalar> EFLU8;
typedef Eigen_Fixed_LU_Decomposition<9, LU_Decomposition_Factory::Scalar> EFLU9;
typedef Eigen_Fixed_LU_Decomposition<10, LU_Decomposition_Factory::Scalar> EFLU10;
/*
typedef Eigen_Fixed_LU_Decomposition<11, LU_Decomposition_Factory::Scalar> EFLU11;
typedef Eigen_Fixed_LU_Decomposition<12, LU_Decomposition_Factory::Scalar> EFLU12;
typedef Eigen_Fixed_LU_Decomposition<13, LU_Decomposition_Factory::Scalar> EFLU13;
typedef Eigen_Fixed_LU_Decomposition<14, LU_Decomposition_Factory::Scalar> EFLU14;
typedef Eigen_Fixed_LU_Decomposition<15, LU_Decomposition_Factory::Scalar> EFLU15;
typedef Eigen_Fixed_LU_Decomposition<16, LU_Decomposition_Factory::Scalar> EFLU16;
typedef Eigen_Fixed_LU_Decomposition<17, LU_Decomposition_Factory::Scalar> EFLU17;
typedef Eigen_Fixed_LU_Decomposition<18, LU_Decomposition_Factory::Scalar> EFLU18;
typedef Eigen_Fixed_LU_Decomposition<19, LU_Decomposition_Factory::Scalar> EFLU19;
typedef Eigen_Fixed_LU_Decomposition<20, LU_Decomposition_Factory::Scalar> EFLU20;
*/

using namespace std;

LU_Decomposition_Factory::
LU_Decomposition_Factory()
{
}

shared_ptr<LU_Decomposition<LU_Decomposition_Factory::Scalar> > LU_Decomposition_Factory::
get_solver(int size,
           Type type) const
{
    switch (type)
    {
    case Type::DEFAULT:
        return get_default(size);
    case Type::EIGEN:
        return make_shared<Eigen_Dense_LU_Decomposition<double> >(size);
    case Type::EIGEN_FIXED:
        return get_fixed_eigen_solver(size);
    case Type::EPETRA:
        return make_shared<Epetra_Dense_LU_Decomposition>(size);
    }
}

shared_ptr<LU_Decomposition<LU_Decomposition_Factory::Scalar> > LU_Decomposition_Factory::
get_default(int size) const
{
    if (size <= 5)
    {
        return get_solver(size,
                          Type::DIRECT);
    }
    else if (size <= 10)
    {
        return get_solver(size,
                          Type::EIGEN_FIXED);
    }
    else if (size < 20)
    {
        return get_solver(size,
                          Type::EIGEN);
    }
    else
    {
        return get_solver(size,
                          Type::EPETRA);
    }
}

shared_ptr<LU_Decomposition<LU_Decomposition_Factory::Scalar> > LU_Decomposition_Factory::
get_fixed_eigen_solver(int size) const
{
    switch (size)
    {
    case 1:
        return make_shared<EFLU1>();
    case 2:
        return make_shared<EFLU2>();
    case 3:
        return make_shared<EFLU3>();
    case 4:
        return make_shared<EFLU4>();
    case 5:
        return make_shared<EFLU5>();
    case 6:
        return make_shared<EFLU6>();
    case 7:
        return make_shared<EFLU7>();
    case 8:
        return make_shared<EFLU8>();
    case 9:
        return make_shared<EFLU9>();
    case 10:
        return make_shared<EFLU10>();
    /*
    case 11:
        return make_shared<EFLU11>();
    case 12:
        return make_shared<EFLU12>();
    case 13:
        return make_shared<EFLU13>();
    case 14:
        return make_shared<EFLU14>();
    case 15:
        return make_shared<EFLU15>();
    case 16:
        return make_shared<EFLU16>();
    case 17:
        return make_shared<EFLU17>();
    case 18:
        return make_shared<EFLU18>();
    case 19:
        return make_shared<EFLU19>();
    case 20:
        return make_shared<EFLU20>();
    */
    default:
        AssertMsg(false, "fixed eigen solver of size (" << size << ") not available");
    }
}
