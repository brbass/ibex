#ifndef Weak_RBF_Sweep_hh
#define Weak_RBF_Sweep_hh

class Weak_RBF_Sweep : public Sweep_Operator
{
    struct Options
    {
        enum class Matrix_Solver
        {
            AMESOS,
            AZTEC,
            EIGEN
        };
    };
};

#endif
