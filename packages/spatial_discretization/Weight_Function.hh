#ifndef Weight_Function_hh
#define Weight_Function_hh

#include "Point.hh"

class Weight_Function : public Point
{
public:
    Weight_Function(int index,
                    int dimension,
                    shared_ptr<Material> material,
                    vector<double> const &position,
                    shared_ptr<RBF_Function> rbf,
                    vector<shared_ptr<Basis_Function> > basis_functions);

private:
    shared_ptr<RBF_Function> rbf_function_;
    vector<shared_ptr<Basis_Function> > basis_functions_;

}

#endif
