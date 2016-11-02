#ifndef Transport_Problem_hh
#define Transport_Problem_hh

#include <memory>
#include <vector>

#include "pugixml.hh"

class Solver;

using std::shared_ptr;
using std::vector;

/*
  High-level class to run a transport problem
*/
class Transport_Problem
{
public:

    // Type of solution
    enum class Problem_Type
    {
        STEADY_STATE,
        K_EIGENVALUE,
        TIME_DEPENDENT
    };

    // Creator
    Transport_Problem(Problem_Type problem_type,
                      shared_ptr<Solver> solver);

    // Solve transport problem
    void solve();

    // Output data to XML file
    void output(pugi::xml_node &output_node) const;

private:
    
    Problem_Type problem_type_;
    shared_ptr<Solver> solver_;

    double k_eigenvalue_;
    vector<double> phi_;
    vector<double> psi_;
};

#endif
