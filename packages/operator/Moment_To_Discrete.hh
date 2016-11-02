#ifndef Moment_To_Discrete_hh
#define Moment_To_Discrete_hh

#include <memory>
#include <vector>

#include "Vector_Operator.hh"

class Angular_Discretization;
class Energy_Discretization;
class Spatial_Discretization;

using std::shared_ptr;
using std::vector;

/*
  Converts moments of the angular flux to the discrete angular flux
*/
class Moment_To_Discrete: public Vector_Operator
{
public:

    // Constructor
    Moment_To_Discrete(shared_ptr<Spatial_Discretization> spatial_discretization,
                       shared_ptr<Angular_Discretization> angular_discretization,
                       shared_ptr<Energy_Discretization> energy_discretization);
    
    virtual void check_class_invariants() const override;

private:

    virtual void apply(vector<double> &x) const override;

    // Output size
    int get_row_size(shared_ptr<Spatial_Discretization> spatial_discretization,
                     shared_ptr<Angular_Discretization> angular_discretization,
                     shared_ptr<Energy_Discretization> energy_discretization);

    // Input size
    int get_column_size(shared_ptr<Spatial_Discretization> spatial_discretization,
                        shared_ptr<Angular_Discretization> angular_discretization,
                        shared_ptr<Energy_Discretization> energy_discretization);
    
    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
};

#endif
