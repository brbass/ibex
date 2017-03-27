#include "Null_Solver.hh"

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Spatial_Discretization.hh"

Null_Solver::
Null_Solver(int solver_print,
            shared_ptr<Spatial_Discretization> spatial_discretization,
            shared_ptr<Angular_Discretization> angular_discretization,
            shared_ptr<Energy_Discretization> energy_discretization):
    Solver(solver_print),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization)
{
}
    
void Null_Solver::
output(pugi::xml_node &output_node) const
{
    pugi::xml_node source = output_node.append_child("krylov_iteration");
    
    spatial_discretization_->output(output_node);
    angular_discretization_->output(output_node);
    energy_discretization_->output(output_node);
}
