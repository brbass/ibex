#include "Preconditioner.hh"

#include <iostream>

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Spatial_Discretization.hh"
#include "Sweep_Operator.hh"

using namespace std;

Preconditioner::
Preconditioner(int row_size,
               int column_size,
               shared_ptr<Spatial_Discretization> spatial_discretization,
               shared_ptr<Angular_Discretization> angular_discretization,
               shared_ptr<Energy_Discretization> energy_discretization,
               shared_ptr<Sweep_Operator> sweeper):
    Vector_Operator(row_size,
                    column_size),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    sweeper_(sweeper)
{
}

