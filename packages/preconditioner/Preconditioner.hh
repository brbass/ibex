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

/*
  Pure virtual class for a preconditioner
*/
class Preconditioner : public Vector_Operator
{
public:
    
    // Constructor
    Preconditioner(int row_size,
                   int column_size,
                   std::shared_ptr<Spatial_Discretization> spatial_discretization,
                   std::shared_ptr<Angular_Discretization> angular_discretization,
                   std::shared_ptr<Energy_Discretization> energy_discretization,
                   std::shared_ptr<Sweep_Operator> sweeper);

    std::shared_ptr<Sweep_Operator> sweeper()
    {
        return sweeper_;
    }
    
protected:

    virtual void apply(std::vector<double> &x) const override = 0;
    
    std::shared_ptr<Spatial_Discretization> spatial_discretization_;
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;
    std::shared_ptr<Sweep_Operator> sweeper_;
};

#endif
