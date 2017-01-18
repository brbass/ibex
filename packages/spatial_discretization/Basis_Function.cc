#include "Basis_Function.hh"

#include "Cartesian_Plane.hh"
#include "Meshless_Function.hh"

using namespace std;

Basis_Function::
Basis_Function(int index,
               int dimension,
               shared_ptr<Meshless_Function> meshless_function,
               vector<shared_ptr<Cartesian_Plane> > boundary_surfaces):
    index_(index),
    dimension_(dimension),
    number_of_boundary_surfaces_(boundary_surfaces.size()),
    radius_(meshless_function->radius()),
    position_(meshless_function->position()),
    meshless_function_(meshless_function),
    boundary_surfaces_(boundary_surfaces)
{
}

