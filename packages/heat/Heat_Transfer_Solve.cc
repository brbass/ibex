#include "Heat_Transfer_Solve.hh"

#include "Amesos.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Vector.h"

#include "Check.hh"
#include "Heat_Transfer_Integration.hh"
#include "Heat_Transfer_Solution.hh"
#include "Weak_Spatial_Discretization.hh"

using namespace std;

Heat_Transfer_Solve::
Heat_Transfer_Solve(shared_ptr<Heat_Transfer_Integration> integration,
                    shared_ptr<Weak_Spatial_Discretization> spatial):
    integration_(integration),
    spatial_(spatial)
{
}

shared_ptr<Heat_Transfer_Solution> Heat_Transfer_Solve::
solve()
{
    // Get needed spatial discretization information
    int number_of_points = spatial_->number_of_points();
    vector<int> const number_of_basis_functions = spatial_->number_of_basis_functions();

    // Get communication classes
    shared_ptr<Epetra_SerialComm> comm
        = make_shared<Epetra_SerialComm>();
    shared_ptr<Epetra_Map> map
        = make_shared<Epetra_Map>(number_of_points, 0, *comm);

    // Get vectors
    shared_ptr<Epetra_Vector> lhs
        = make_shared<Epetra_Vector>(*map);
    shared_ptr<Epetra_Vector> rhs
        = make_shared<Epetra_Vector>(Copy,
                                     *map,
                                     &integration_->rhs()[0]);
    
    // Get matrix
    shared_ptr<Epetra_CrsMatrix> mat
        = make_shared<Epetra_CrsMatrix>(Copy,
                                        *map,
                                        &number_of_basis_functions[0],
                                        true);
    for (int i = 0; i < number_of_points; ++i)
    {
        vector<int> const indices = spatial_->weight(i)->basis_function_indices();
        vector<double> const values = integration_->matrix()[i];
        Assert(indices.size() == number_of_basis_functions[i]);
        Assert(values.size() == number_of_basis_functions[i]);
        
        mat->InsertGlobalValues(i,
                                number_of_basis_functions[i],
                                &values[0],
                                &indices[0]);
    }
    mat->FillComplete();
    mat->OptimizeStorage();
    
    // Get problem and solver
    shared_ptr<Epetra_LinearProblem> problem
        = make_shared<Epetra_LinearProblem>(mat.get(),
                                            lhs.get(),
                                            rhs.get());
    Amesos factory;
    shared_ptr<Amesos_BaseSolver> solver(factory.Create("Klu",
                                                        *problem));
    AssertMsg(solver->SymbolicFactorization() == 0, "Amesos solver symbolic factorization failed");
    AssertMsg(solver->NumericFactorization() == 0, "Amesos solver numeric factorization failed");

    // Solve problem
    AssertMsg(solver->Solve() == 0, "Heat solver failed");

    // Get solution
    vector<double> coefficients(number_of_points);
    for (int i = 0; i < number_of_points; ++i)
    {
        coefficients[i] = (*lhs)[i];
        cout << coefficients[i] << endl;
    }

    return make_shared<Heat_Transfer_Solution>(spatial_,
                                               coefficients);
}
