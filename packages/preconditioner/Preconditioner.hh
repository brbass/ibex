#ifndef Preconditioner_hh
#define Preconditioner_hh

#include <memory>
#include <vector>

#include "Vector_Operator.hh"

#include "Sweep_Operator.hh"

class Angular_Discretization;
class Energy_Discretization;
class Nuclear_Data;
class Source_Data;
class Spatial_Discretization;

using std::shared_ptr;
using std::vector;

/*
  Pure virtual class for a preconditioner
*/
class Preconditioner : public Vector_Operator
{
public:
    
    // Constructor
    Preconditioner(int row_size,
                   int column_size,
                   shared_ptr<Spatial_Discretization> spatial_discretization,
                   shared_ptr<Angular_Discretization> angular_discretization,
                   shared_ptr<Energy_Discretization> energy_discretization,
                   shared_ptr<Sweep_Operator> sweeper);

    shared_ptr<Sweep_Operator> sweeper()
    {
        return sweeper_;
    }
    
protected:

    virtual void apply(vector<double> &x) const override = 0;
    
    shared_ptr<Spatial_Discretization> spatial_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
    shared_ptr<Sweep_Operator> sweeper_;
};

#endif
