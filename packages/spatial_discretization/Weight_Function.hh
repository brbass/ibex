#ifndef Weight_Function_hh
#define Weight_Function_hh

#include "Point.hh"

class Basis_Function;
class Cartesian_Plane;
class RBF_Function;

class Weight_Function : public Point
{
public:
    
    Weight_Function(int index,
                    int dimension,
                    vector<double> const &position,
                    shared_ptr<Meshless_Function> meshless_function,
                    vector<shared_ptr<Basis_Function> > basis_functions,
                    shared_ptr<Solid_Geometry> solid_geometry,
                    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces);

private:
    
    int index_;
    int dimension_;
    vector<double> position_;
    shared_ptr<Meshless_Function> meshless_function_;
    vector<shared_ptr<Basis_Function> > basis_functions_;
    shared_ptr<Solid_Geometry> solid_geometry_;
    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces_;

}

#endif
