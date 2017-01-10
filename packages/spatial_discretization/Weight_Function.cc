#include "Weight_Function.hh"

#include "Basis_Function.hh"
#include "Cartesian_Plane.hh"
#include "Meshless_Function.hh"
#include "Solid_Geometry.hh"

Weight_Function::
Weight_Function(int index,
                int dimension,
                vector<double> const &position,
                shared_ptr<Meshless_Function> meshless_function,
                vector<shared_ptr<Basis_Function> > basis_functions,
                shared_ptr<Solid_Geometry> solid_geometry,
                 vector<shared_ptr<Cartesian_Plane> > boundary_surfaces):
    index_(index),
    dimension_(dimension),
    solid_geometry_(solid_geometry),
    position_(position),
    meshless_function_(meshless_function),
    basis_functions_(basis_functions),
    boundary_surfaces_(boundary_surfaces)
{
    
}
