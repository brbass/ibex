#include <memory>
#include <mpi.h>
#include <vector>

#include "Energy_Discretization.hh"
#include "LDFE_Quadrature.hh"
#include "VERA_Solid_Geometry.hh"

using namespace std;

int main(int argc, char **argv)
{
    // Initialize MPI
    MPI_Init(&argc, &argv);

    shared_ptr<Energy_Discretization> energy
        = make_shared<Energy_Discretization>(2);
    shared_ptr<Angular_Discretization> angular
        = make_shared<LDFE_Quadrature>(2,
                                       2,
                                       1);
    shared_ptr<VERA_Solid_Geometry> solid
        = make_shared<VERA_Solid_Geometry>(false,
                                           [](vector<double> const &){return 600;},
                                           angular,
                                           energy);
        
    // Close MPI
    MPI_Finalize();
}
