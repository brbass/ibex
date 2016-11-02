#ifndef Dense_Solve_hh
#define Dense_Solve_hh

#include <memory>
#include <vector>

class Trilinos_Dense_Solve;

/*
  Generalized class for "dumb" dense matrix solution
*/
class Dense_Solve
{

public:

    // Creator
    Dense_Solve(unsigned size);

    // Solve problem Ax=b
    void solve(std::vector<double> &a_data, 
               std::vector<double> &b_data,
               std::vector<double> &x_data);

    // Size of square matrix
    unsigned size()
    {
        return size_;
    }

private:
    
    unsigned size_;
    
    std::shared_ptr<Trilinos_Dense_Solve> trilinos_solver_;
};

#endif
