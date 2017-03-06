#include "RBF_Collocation_Sweep.hh"

#include <cmath>
#include <limits>
#include <memory>
#include <vector>

#include <Amesos.h>
#include <AztecOO.h>
#include <AztecOO_ConditionNumber.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseSolver.h>
#include <Epetra_SerialDenseVector.h>
#include <Ifpack.h>
#include <Ifpack_Preconditioner.h>
// #include <mpi.h>

#include "Angular_Discretization.hh"
#include "Boundary_Source.hh"
#include "Check_Equality.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "RBF_Discretization.hh"
#include "Spatial_Discretization.hh"
#include "Symmetric_Sparse_Storage.hh"
#include "Timer.hh"
#include "Transport_Discretization.hh"
#include "XML_Functions.hh"

using namespace std;
namespace ce = Check_Equality;

RBF_Collocation_Sweep::
RBF_Collocation_Sweep(Solution_Variable solution_variable,
                      Matrix_Solver matrix_solver,
                      Condition_Calculation condition_calculation,
                      shared_ptr<RBF_Discretization> rbf_discretization,
                      shared_ptr<Angular_Discretization> angular_discretization,
                      shared_ptr<Energy_Discretization> energy_discretization,
                      shared_ptr<Transport_Discretization> transport_discretization):
    Sweep_Operator(Sweep_Type::ORDINATE,
                   transport_discretization),
    solution_variable_(solution_variable),
    matrix_solver_(matrix_solver),
    condition_calculation_(condition_calculation),
    rbf_discretization_(rbf_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization)
{
    int number_of_ordinates = angular_discretization->number_of_ordinates();
    int number_of_groups = energy_discretization->number_of_groups();
    reflection_tolerance_ = 1000 * numeric_limits<double>::epsilon();
    num_calls_.assign(number_of_ordinates * number_of_groups, 0);
    num_iterations_.assign(number_of_ordinates * number_of_groups, 0);
    num_averages_.assign(number_of_ordinates * number_of_groups, 0);
    average_local_condition_numbers_.assign(number_of_ordinates * number_of_groups, 0);
    
    initialize_trilinos();
    calculate_condition_numbers();
    check_class_invariants();
}

void RBF_Collocation_Sweep::
apply(vector<double> &x) const
{
    int dimension = rbf_discretization_->dimension();
    int number_of_points = rbf_discretization_->number_of_points();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_boundary_points = rbf_discretization_->number_of_boundary_points();
    int number_of_internal_points = rbf_discretization_->number_of_internal_points();
    int number_of_neighbors = rbf_discretization_->number_of_neighbors();
    int number_of_augments = transport_discretization_->number_of_augments();
    int psi_size = transport_discretization_->psi_size();
    vector<int> const boundary_points = rbf_discretization_->boundary_points();
    vector<int> const internal_points = rbf_discretization_->internal_points();
    
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        vector<double> const direction = angular_discretization_->direction(o);
        
        for (int g = 0; g < number_of_groups; ++g)
        {
            for (int b = 0; b < number_of_boundary_points; ++b)
            {
                int i = boundary_points[b];
                
                double sum = 0;

                vector<double> const boundary_normal
                    = rbf_discretization_->rbf_point(i)->normal();
                
                for (int d = 0; d < dimension; ++d)
                {
                    int k_bn = d + dimension * b;
                    int k_ord = d + dimension * o;
                    
                    sum += boundary_normal[d] * direction[d];
                }
                
                if (sum < 0)
                {
                    set_boundary_rhs(b, i, o, g, x);
                }
                else
                {
                    set_internal_rhs(i, o, g, x);
                }
            }
            
            for (int p = 0; p < number_of_internal_points; ++p)
            {
                int i = internal_points[p];
                
                set_internal_rhs(i, o, g, x);
            }
            
            // Perform matrix solve

            int k = g + number_of_groups * o;

            num_calls_[k] += 1;
            switch (matrix_solver_)
            {
            case Matrix_Solver::AMESOS:
                (*amesos_solver_[k])->Solve();
                break;
            case Matrix_Solver::AZTEC:
                aztec_solver_[k]->Iterate(max_iterations_, tolerance_);
                num_iterations_[k] += aztec_solver_[k]->NumIters();
                break;
            }
            
            // Update values for this group and ordinate
            
            for (int i = 0; i < number_of_points; ++i)
            {
                vector<int> const neighbors
                    = rbf_discretization_->rbf_point(i)->neighbor_indices();

                int k_x = g + number_of_groups * (o + number_of_ordinates * i);

                switch(solution_variable_)
                {
                case Solution_Variable::COEFFICIENT:
                {
                    double value = 0;
                    for (int n = 0; n < number_of_neighbors; ++n)
                    {
                        int j = neighbors[n];
                        
                        value += (*lhs_)[j] * rbf_discretization_->basis(i,
                                                                         j,
                                                                         g,
                                                                         o);
                    }
                    
                    x[k_x] = value;
                    break;
                }
                case Solution_Variable::PSI:
                {
                    x[k_x] = (*lhs_)[i];
                    break;
                }
                }
            }
        }
    }
    
    // Update augments
    
    for (int b = 0; b < number_of_boundary_points; ++b)
    {
        int i = boundary_points[b];
        
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k_b = psi_size + g + number_of_groups * (o + number_of_ordinates * b);
                int k_psi = g + number_of_groups * (o + number_of_ordinates * i);
                    
                x[k_b] = x[k_psi];
            }
        }
    }
}

void RBF_Collocation_Sweep::
initialize_trilinos()
{
    int number_of_points = rbf_discretization_->number_of_points();
    int number_of_neighbors = rbf_discretization_->number_of_neighbors();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    comm_ = make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    map_ = make_shared<Epetra_Map>(number_of_points, 0, *comm_);
    lhs_ = make_shared<Epetra_Vector>(*map_);
    rhs_ = make_shared<Epetra_Vector>(*map_);
    
    lhs_->PutScalar(1.0);
    rhs_->PutScalar(1.0);
    
    mat_.resize(number_of_groups * number_of_ordinates);
    problem_.resize(number_of_groups * number_of_ordinates);

    // Fill the entries of the matrices with ones and optimize storage
    vector<double> ones(number_of_neighbors, 1);
    for (int g = 0; g < number_of_groups; ++g)
    {
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            int k = g + number_of_groups * o;
            
            mat_[k] = make_shared<Epetra_CrsMatrix>(Copy, *map_, number_of_neighbors, true);
            
            for (int i = 0; i < number_of_points; ++i)
            {
                vector<int> const neighbors = rbf_discretization_->rbf_point(i)->neighbor_indices();
                mat_[k]->InsertGlobalValues(i, number_of_neighbors, &ones[0], &neighbors[0]);
            }
            mat_[k]->FillComplete();
            mat_[k]->OptimizeStorage();
            
            problem_[k] = make_shared<Epetra_LinearProblem>(mat_[k].get(),
                                                            lhs_.get(),
                                                            rhs_.get());
        }
    }

    // Fill the matrices with the correct values
    switch(solution_variable_)
    {
    case Solution_Variable::PSI:
    {
        for (int i = 0; i < number_of_points; ++i)
        {
            shared_ptr<RBF_Function> rbf_function = rbf_discretization_->rbf_point(i)->rbf_function();
            
            if (rbf_function->direction_dependent())
            {
                if (rbf_function->energy_dependent())
                {
                    for (int o = 0; o < number_of_ordinates; ++o)
                    {
                        for (int g = 0; g < number_of_groups; ++g)
                        {
                            update_local_matrix(i, o, g);
                            check_local_matrix(i, o, g);
                            set_point(i, o, g);
                        }
                    }                    
                }
                else
                {
                    for (int o = 0; o < number_of_ordinates; ++o)
                    {
                        update_local_matrix(i, o, 0);
                        check_local_matrix(i, o, 0);
                        for (int g = 0; g < number_of_groups; ++g)
                        {
                            set_point(i, o, g);
                        }
                    }                    
                }
            }
            else if (rbf_function->energy_dependent())
            {
                for (int g = 0; g < number_of_groups; ++g)
                {
                    update_local_matrix(i, 0, g);
                    check_local_matrix(i, 0, g);
                    for (int o = 0; o < number_of_ordinates; ++o)
                    {
                        set_point(i, o, g);
                    }
                }
            }
            else
            {
                update_local_matrix(i, 0, 0);  
                check_local_matrix(i, 0, 0);
                for (int g = 0; g < number_of_groups; ++g)
                {
                    for (int o = 0; o < number_of_ordinates; ++o)
                    {
                        set_point(i, o, g);
                    }
                }
            }
        }
        break;
    }
    case Solution_Variable::COEFFICIENT:
    {
        for (int i = 0; i < number_of_points; ++i)
        {
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                for (int g = 0; g < number_of_groups; ++g)
                {
                    set_point(i, o, g);
                }
            }
        }
        break;
    }
    }

    switch (matrix_solver_)
    {
    case Matrix_Solver::AMESOS:
        amesos_solver_.resize(number_of_groups * number_of_ordinates);
        
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k = g + number_of_groups * o;
                
                Amesos factory;
                // Serial: Klu, Lapack, Umfpack
                // Parallel: Mumps, Superludist
                amesos_solver_[k] = make_shared<Amesos_BaseSolver*>(factory.Create("Klu", *problem_[k]));
                
                if (*amesos_solver_[k] == NULL)
                {
                    AssertMsg(false, "specified solver is not available");
                }
                
                (*amesos_solver_[k])->SymbolicFactorization();
                (*amesos_solver_[k])->NumericFactorization();
            }
        }
        
        break;
    case Matrix_Solver::AZTEC:
        max_iterations_ = 1000;
        tolerance_ = 1e-6;
        num_iterations_.assign(number_of_groups * number_of_ordinates, 0);
        aztec_solver_.resize(number_of_groups * number_of_ordinates);
        
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k = g + number_of_groups * o;

                aztec_solver_[k] = make_shared<AztecOO>(*problem_[k]);

                aztec_solver_[k]->SetAztecOption(AZ_solver, AZ_gmres);
                aztec_solver_[k]->SetAztecOption(AZ_kspace, 100);
                aztec_solver_[k]->SetAztecOption(AZ_precond, AZ_none);
                aztec_solver_[k]->SetAztecOption(AZ_output, AZ_none);
            }
        }
        
        break;
    }
}

void RBF_Collocation_Sweep::
set_point(int i,
          int o,
          int g)
{
    shared_ptr<RBF_Point> rbf_point = rbf_discretization_->rbf_point(i);

    switch(rbf_point->point_type())
    {
    case RBF_Point::Point_Type::INTERNAL:
    {
        set_internal_point(i, o, g);
        break;
    }
    case RBF_Point::Point_Type::BOUNDARY:
    {
        int dimension = rbf_discretization_->dimension();
        vector<double> const boundary_normal = rbf_point->normal();
        vector<double> const ordinates = angular_discretization_->ordinates();
        vector<double> const direction = angular_discretization_->direction(o);
        
        double sum = 0;
        for (int d = 0; d < dimension; ++d)
        {
            sum += boundary_normal[d] * direction[d];
        }
        
        if (sum < 0)
        {
            set_boundary_point(i, o, g);
        }
        else
        {
            set_internal_point(i, o, g);
        }
        break;
    }
    }
}

void RBF_Collocation_Sweep::
set_boundary_point(int i,
                   int o,
                   int g)
{
    // int dimension = rbf_discretization_->dimension();
    // int number_of_points = rbf_discretization_->number_of_points();
    int number_of_groups = energy_discretization_->number_of_groups();
    // int number_of_ordinates = angular_discretization_->number_of_ordinates();
    // int number_of_augments = transport_discretization_->number_of_augments();
    int number_of_neighbors = rbf_discretization_->number_of_neighbors();
    
    vector<int> const neighbors = rbf_discretization_->rbf_point(i)->neighbor_indices();
    
    // Replace matrix values
    vector<double> data(number_of_neighbors, 0);

    switch(solution_variable_)
    {
    case Solution_Variable::COEFFICIENT:
    {
        for (int n = 0; n < number_of_neighbors; ++n)
        {
            int j = neighbors[n];
            
            data[n] = rbf_discretization_->basis(i,
                                                 j,
                                                 g,
                                                 o);
        }
        break;
    }
    case Solution_Variable::PSI:
    {
        Check(neighbors[0] == i);
        
        data[0] = 1;
        break;
    }
    }
    
    int k = g + number_of_groups * o;
    
    mat_[k]->ReplaceGlobalValues(i,
                                 number_of_neighbors,
                                 &data[0],
                                 &neighbors[0]);
}

void RBF_Collocation_Sweep::
set_boundary_rhs(int b,
                 int i,
                 int o,
                 int g,
                 vector<double> const &x) const
{
    // int dimension = rbf_discretization_->dimension();
    // int number_of_points = rbf_discretization_->number_of_points();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    // int number_of_augments = transport_discretization_->number_of_augments();
    int psi_size = transport_discretization_->psi_size();

    shared_ptr<RBF_Point> point
        = rbf_discretization_->rbf_point(i);
    shared_ptr<Boundary_Source> source
        = point->boundary_source();
    vector<double> const boundary_source = source->boundary_source();
    vector<double> const alpha = source->alpha();
    
    vector<double> const boundary_normal = point->normal();
    
    double rhs = 0;
    if (alpha[g] > reflection_tolerance_)
    {
        int o1 = angular_discretization_->reflect_ordinate(o,
                                                           boundary_normal);
        int k_ref = psi_size + g + number_of_groups * (o1 + number_of_ordinates * b);
        
        rhs = alpha[b] * x[k_ref];
    }
    
    if (include_boundary_source_)
    {
        int k_bs = g + number_of_groups * o;
        
        rhs += boundary_source[k_bs];
    }
    
    (*rhs_)[i] = rhs;
}

void RBF_Collocation_Sweep::
set_internal_point(int i,
                   int o,
                   int g)
{
    int dimension = rbf_discretization_->dimension();
    // int number_of_points = rbf_discretization_->number_of_points();
    int number_of_groups = energy_discretization_->number_of_groups();
    // int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_neighbors = rbf_discretization_->number_of_neighbors();
    // int number_of_augments = transport_discretization_->number_of_augments();
    // int psi_size = transport_discretization_->psi_size();
    
    shared_ptr<RBF_Point> rbf_point = rbf_discretization_->rbf_point(i);
    vector<int> const neighbors = rbf_point->neighbor_indices();
    vector<double> const sigma_t = rbf_point->material()->sigma_t()->data();
    vector<double> const direction = angular_discretization_->direction(o);
    
    // Replace matrix values
    
    vector<double> data(number_of_neighbors, 0);
    
    for (int n = 0; n < number_of_neighbors; ++n)
    {
        int j = neighbors[n];
        
        vector<double> const gradient
            = rbf_discretization_->gradient_basis(i,
                                                  j,
                                                  g,
                                                  o);
        
        double derivative = 0;
        for (int d = 0; d < dimension; ++d)
        {
            derivative += direction[d] * gradient[d];
        }
        
        double basis = rbf_discretization_->basis(i,
                                                  j,
                                                  g,
                                                  o);
        
        data[n] = derivative + sigma_t[g] * basis;
    }
    
    switch(solution_variable_)
    {
    case Solution_Variable::COEFFICIENT:
        break;
    case Solution_Variable::PSI:
        convert_to_psi(i,
                       g,
                       o,
                       data);
        break;
    }
    
    int k = g + number_of_groups * o;
    
    mat_[k]->ReplaceGlobalValues(i,
                                 number_of_neighbors,
                                 &data[0],
                                 &neighbors[0]);
}

void RBF_Collocation_Sweep::
check_local_matrix(int i,
                   int o, 
                   int g)
{
    if (false)
    {
        switch(solution_variable_)
        {
        case Solution_Variable::COEFFICIENT:
            break;
        case Solution_Variable::PSI:
            int dimension = rbf_discretization_->dimension();
            int number_of_points = rbf_discretization_->number_of_points();
            int number_of_groups = energy_discretization_->number_of_groups();
            int number_of_ordinates = angular_discretization_->number_of_ordinates();
            int number_of_neighbors = rbf_discretization_->number_of_neighbors();
            int number_of_augments = transport_discretization_->number_of_augments();
            int psi_size = transport_discretization_->psi_size();
    
            shared_ptr<RBF_Point> rbf_point = rbf_discretization_->rbf_point(i);
            vector<int> const neighbors = rbf_point->neighbor_indices();
            vector<double> const sigma_t = rbf_point->material()->sigma_t();
            vector<double> const equation_position = rbf_point->position();
            vector<double> const direction = angular_discretization_->direction(o);
    
            // Get test values
    
            vector<double> data(number_of_neighbors, 0);
    
            for (int n = 0; n < number_of_neighbors; ++n)
            {
                int j = neighbors[n];
        
                shared_ptr<RBF_Point> basis_point = rbf_discretization_->rbf_point(j);
                shared_ptr<RBF_Function> basis_rbf = rbf_point->rbf_function();
                vector<double> const shape_parameter = basis_point->shape_parameter();
                vector<double> const basis_position = basis_point->position();
                double dist = rbf_discretization_->distances(g)->get_element(i, j);
        
                double basis = basis_rbf->basis(g,
                                                shape_parameter[g],
                                                dist,
                                                equation_position,
                                                basis_position,
                                                direction);
        
                data[n] = basis;
            }
    
            convert_to_psi(i,
                           g,
                           o,
                           data);
    
            vector<double> check_data(number_of_neighbors, 0);
            check_data[0] = 1.;
            
            if (!ce::approx(check_data, data, 1e-10))
            {
                cout << "Test values do not match for local RBF matrix:\t";
                cout << check_data[0] << " / " << data[0] << "  ";
                cout << check_data[1] << " / " << data[1] << endl;
            }

            break;
        }
    }
}

void RBF_Collocation_Sweep::
set_internal_rhs(int i,
                 int o,
                 int g,
                 vector<double> const &x) const
{
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    
    int k_x = g + number_of_groups * (o + number_of_ordinates * i);
    
    (*rhs_)[i] = x[k_x];
}

void RBF_Collocation_Sweep::
update_local_matrix(int i,
                    int o,
                    int g)
{
    // Create matrix transposed to convert coefficients
    // Usually, basis function would depend on column
    // Because of transpose, basis function depends on row
    // Could alternately solve using transpose. 

    int number_of_neighbors = rbf_discretization_->number_of_neighbors();
    vector<int> const neighbors = rbf_discretization_->rbf_point(i)->neighbor_indices();
    
    // Create matrix and solver if they don't exist
    if (!local_mat_)
    {
        local_mat_ = make_shared<Epetra_SerialDenseMatrix>(number_of_neighbors,
                                                           number_of_neighbors);
    }
    if (!local_solver_)
    {
        local_solver_ = make_shared<Epetra_SerialDenseSolver>();
    }
    
    // Fill matrix
    for (int n1 = 0; n1 < number_of_neighbors; ++n1)
    {
        int j = neighbors[n1];
        
        for (int n2 = 0; n2 < number_of_neighbors; ++n2)
        {
            int k = neighbors[n2];
            
            (*local_mat_)(n1, n2) = rbf_discretization_->basis(k,
                                                               j,
                                                               g,
                                                               o);
        }
    }

    local_solver_->SetMatrix(*local_mat_);
    local_solver_->Factor();
}

void RBF_Collocation_Sweep::
convert_to_psi(int i,
               int g,
               int o,
               vector<double> &data) const
{
    int number_of_neighbors = rbf_discretization_->number_of_neighbors();
    
    Epetra_SerialDenseVector x(number_of_neighbors);
    Epetra_SerialDenseVector b(View, &data[0], number_of_neighbors);
    
    local_solver_->SetVectors(x, b);
    local_solver_->Solve();

    switch (condition_calculation_)
    {
    case Condition_Calculation::NONE:
        break;
    default:
    {
        double condition;
        local_solver_->ReciprocalConditionEstimate(condition);
        
        int k = g + energy_discretization_->number_of_groups() * o;
        average_local_condition_numbers_[k] = (average_local_condition_numbers_[k] + 1 / condition) / (num_averages_[k] + 1);
        num_averages_[k] += 1;
        break;
    }
    }
}

void RBF_Collocation_Sweep::
calculate_condition_numbers()
{
    int max_iterations = 1000;
    double tolerance = 1e-10;
    
    condition_numbers_.resize(problem_.size());
    AztecOOConditionNumber aztec_condition;
    Ifpack ifpack;
    shared_ptr<Ifpack_Preconditioner*> preconditioner;
    for (int i = 0; i < problem_.size(); ++i)
    {
        switch (condition_calculation_)
        {
        case Condition_Calculation::NONE:
            break;
        case Condition_Calculation::CHEAP:
            preconditioner = make_shared<Ifpack_Preconditioner*>(ifpack.Create("Amesos", mat_[i].get()));
            (*preconditioner)->Initialize();
            (*preconditioner)->Compute();
            condition_numbers_[i] = (*preconditioner)->Condest(Ifpack_Cheap);
            break;
        case Condition_Calculation::EXPENSIVE:
            preconditioner = make_shared<Ifpack_Preconditioner*>(ifpack.Create("Amesos", mat_[i].get()));
            (*preconditioner)->Initialize();
            (*preconditioner)->Compute();
            condition_numbers_[i] = (*preconditioner)->Condest(Ifpack_GMRES);
            break;
        case Condition_Calculation::AZTEC:
        {
            aztec_condition.initialize(*mat_[i]);
            int status = aztec_condition.computeConditionNumber(2000, 1e-6);
            if (status != 0)
            {
                cout << "status from AztecOO computeConditionNumber: " << status << endl;
            }
            condition_numbers_[i] = aztec_condition.getConditionNumber();
        }
        }
    }
}

void RBF_Collocation_Sweep::
output(pugi::xml_node output_node) const
{
    pugi::xml_node sweep = output_node.append_child("sweep_operator");
    
    XML_Functions::append_child(sweep, "rbf_collocation", "sweep_type");

    switch (condition_calculation_)
    {
    case Condition_Calculation::NONE:
        break;
    default:
        XML_Functions::append_child(sweep, condition_numbers_, "condition_numbers", "group-ordinate");

        switch (solution_variable_)
        {
        case Solution_Variable::PSI:
            XML_Functions::append_child(sweep, average_local_condition_numbers_, "local_average_condition_numbers", "group-ordinate");
        case Solution_Variable::COEFFICIENT:
            break;
        }
        break;
    }
    
    switch (solution_variable_)
    {
    case Solution_Variable::PSI:
        XML_Functions::append_child(sweep, "psi", "solution_variable");
        break;
    case Solution_Variable::COEFFICIENT:
        XML_Functions::append_child(sweep, "coefficient", "solution_variable");
        break;
    }

    switch (condition_calculation_)
    {
    case Condition_Calculation::NONE:
        XML_Functions::append_child(sweep, "none", "condition_calculation");
        break;
    case Condition_Calculation::CHEAP:
        XML_Functions::append_child(sweep, "cheap", "condition_calculation");
        break;
    case Condition_Calculation::EXPENSIVE:
        XML_Functions::append_child(sweep, "expensive", "condition_calculation");
        break;
    case Condition_Calculation::AZTEC:
        XML_Functions::append_child(sweep, "aztec", "condition_calculation");
        break;
    }

    XML_Functions::append_child(sweep, num_calls_, "number_of_matrix_solution_calls", "group-ordinate");

    for (int i = 0; i < num_calls_.size(); ++i)
    {
        num_iterations_[i] = round(1.0 * num_iterations_[i] / num_calls_[i]);
    }
    
    
    switch (matrix_solver_)
    {
    case Matrix_Solver::AZTEC:
        XML_Functions::append_child(sweep, "aztec", "matrix_solver");
        XML_Functions::append_child(sweep, num_iterations_, "average_number_of_iterations");
        break;
    case Matrix_Solver::AMESOS:
        XML_Functions::append_child(sweep, "amesos", "matrix_solver");
        break;
    }
}

void RBF_Collocation_Sweep::
check_class_invariants() const
{
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    Assert(mat_.size() == number_of_groups * number_of_ordinates);
    Assert(problem_.size() == number_of_groups * number_of_ordinates);
    Assert(condition_numbers_.size() == number_of_groups * number_of_ordinates);
    Assert(num_averages_.size() == number_of_groups * number_of_ordinates);
    Assert(average_local_condition_numbers_.size() == number_of_groups * number_of_ordinates);
}
