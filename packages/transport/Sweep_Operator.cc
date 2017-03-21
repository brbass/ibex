#include "Sweep_Operator.hh"

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Spatial_Discretization.hh"
#include "Transport_Discretization.hh"

using std::shared_ptr;

Sweep_Operator::
Sweep_Operator(Sweep_Type sweep_type,
               shared_ptr<Transport_Discretization> transport_discretization):
    Square_Vector_Operator(),
    include_boundary_source_(false),
    sweep_type_(sweep_type),
    transport_discretization_(transport_discretization)
{
    switch(sweep_type)
    {
    case Sweep_Operator::Sweep_Type::MOMENT:
        size_ = (transport_discretization->phi_size()
                 + transport_discretization->number_of_augments());
    case Sweep_Operator::Sweep_Type::ORDINATE:
        size_ = (transport_discretization->psi_size()
                 + transport_discretization->number_of_augments());
    }
}

