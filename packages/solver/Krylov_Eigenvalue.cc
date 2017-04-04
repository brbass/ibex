#include "Krylov_Eigenvalue.hh"

#include <AztecOO.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>

#include "Angular_Discretization.hh"
#include "Convergence_Measure.hh"
#include "Energy_Discretization.hh"
#include "Epetra_Operator_Interface.hh"
#include "Spatial_Discretization.hh"
#include "Transport_Discretization.hh"
#include "Vector_Operator.hh"
#include "XML_Node.hh"

using namespace std;

Krylov_Eigenvalue ::
Krylov_Eigenvalue(Options options,
                    shared_ptr<Spatial_Discretization> spatial_discretization,
                    shared_ptr<Angular_Discretization> angular_discretization,
                    shared_ptr<Energy_Discretization> energy_discretization,
                    shared_ptr<Transport_Discretization> transport_discretization,
                  shared_ptr<Vector_Operator> source_operator,
                  shared_ptr<Vector_Operator> flux_operator,
                  vector<shared_ptr<Vector_Operator> > value_operators):
    Solver(options.solver_print,
           Solver::Type::K_EIGENVALUE),
    options_(options),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    transport_discretization_(transport_discretization),
    source_operator_(source_operator),
    flux_operator_(flux_operator),
    value_operators_(value_operators)
{
}

void Krylov_Eigenvalue::
solve()
{
    
}

void Krylov_Eigenvalue::
output(XML_Node output_node) const
{
    
}

void Krylov_Eigenvalue::
check_class_invariants() const
{
    
}
