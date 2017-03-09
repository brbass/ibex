#include "Weak_RBF_Sweep.hh"

#include <iostream>

#include <Amesos.h>
#include <AztecOO.h>
#include <AztecOO_ConditionNumber.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>
// #include <Epetra_SerialDenseMatrix.h>
// #include <Epetra_SerialDenseSolver.h>
// #include <Epetra_SerialDenseVector.h>
// #include <Ifpack.h>
// #include <Ifpack_Preconditioner.h>

#include "Angular_Discretization.hh"
#include "Basis_Function.hh"
#include "Boundary_Source.hh"
#include "Cartesian_Plane.hh"
#include "Check.hh"
#include "Cross_Section.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "Transport_Discretization.hh"
#include "XML_Node.hh"

using std::make_shared;
using std::shared_ptr;
using std::vector;

Weak_RBF_Sweep::
Weak_RBF_Sweep(Options options,
               shared_ptr<Weak_Spatial_Discretization> spatial_discretization,
               shared_ptr<Angular_Discretization> angular_discretization,
               shared_ptr<Energy_Discretization> energy_discretization,
               shared_ptr<Transport_Discretization> transport_discretization):
    Sweep_Operator(Sweep_Type::ORDINATE,
                   transport_discretization),
    options_(options),
    spatial_discretization_(spatial_discretization),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization)
{
    switch (options_.solver)
    {
    case Options::Solver::AMESOS:
        solver_ = make_shared<Amesos_Solver>(*this);
        break;
    case Options::Solver::AZTEC:
        Aztec_Solver::Options aztec_options;
        solver_ = make_shared<Aztec_Solver>(*this,
                                            aztec_options);
    }

    check_class_invariants();
}

void Weak_RBF_Sweep::
check_class_invariants() const
{
    Assert(spatial_discretization_);
    Assert(angular_discretization_);
    Assert(energy_discretization_);
    Assert(solver_);
}

void Weak_RBF_Sweep::
apply(vector<double> &x) const
{
    solver_->solve(x);
    update_augments(x);
}

void Weak_RBF_Sweep::
update_augments(vector<double> &x) const
{
    int number_of_boundary_points = spatial_discretization_->number_of_boundary_points();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_groups = energy_discretization_->number_of_groups();
    int psi_size = transport_discretization_->psi_size();
    bool has_reflection = transport_discretization_->has_reflection();

    // No augments if there is no reflection
    if (!has_reflection)
    {
        return;
    }

    // Iterate over boundary indices
    for (int b = 0; b < number_of_boundary_points; ++b)
    {
        // Get global index
        int i = spatial_discretization_->boundary_basis(b)->index();

        for (int o = 0; o < number_of_ordinates; ++o)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                // Set boundary augment to current value of angular flux
                int k_b = psi_size + g + number_of_groups * (o + number_of_ordinates * b);
                int k_psi = g + number_of_groups * (o + number_of_ordinates * i);
                x[k_b] = x[k_psi];
            }
        }
    }
}

void Weak_RBF_Sweep::
get_rhs(int i,
        int o,
        int g,
        vector<double> const &x,
        double &value) const
{
    // Get data
    shared_ptr<Weight_Function> weight = spatial_discretization_->weight(i);
    vector<double> const is_b_w = weight->is_b_w();
    vector<double> const direction = angular_discretization_->direction(o);
    int number_of_basis_functions = weight->number_of_basis_functions();
    int number_of_boundary_surfaces = weight->number_of_boundary_surfaces();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_groups = energy_discretization_->number_of_groups();
    int dimension = spatial_discretization_->dimension();
    int psi_size = transport_discretization_->psi_size();
    bool has_reflection = transport_discretization_->has_reflection();

    value = 0;
    // Add reflection and boundary source contribution
    {
        // Get sum of normals and integrals
        vector<double> sum(dimension, 0);
        for (int s = 0; s < number_of_boundary_surfaces; ++s)
        {
            shared_ptr<Cartesian_Plane> surface = weight->boundary_surface(s);
            int surface_dimension = surface->surface_dimension();
            double const normal = surface->normal();
            
            // Only for incoming surfaces
            double dot = normal * direction[surface_dimension];
            if (dot < 0)
            {
                shared_ptr<Boundary_Source> source = surface->boundary_source();
                double local_sum = 0;
                // Add reflection
                if (has_reflection)
                {
                    double const alpha = source->alpha()[g];
                    vector<double> normal_vec(dimension, 0);
                    normal_vec[surface_dimension] = normal;
                    int o_ref = angular_discretization_->reflect_ordinate(o,
                                                                          normal_vec);
                    
                    // Sum contributions of basis functions
                    for (int j = 0; j < number_of_basis_functions; ++j)
                    {
                        shared_ptr<Basis_Function> basis = weight->basis_function(j);
                        if (basis->point_type() == Basis_Function::Point_Type::BOUNDARY)
                        {
                            int aug_index = basis->boundary_index();
                            int is_index = s + number_of_boundary_surfaces * j;
                            int psi_index = psi_size + g + number_of_groups * (o_ref + number_of_ordinates * aug_index);
                            local_sum += is_b_w[is_index] * x[psi_index] * alpha;
                        }
                    }
                }
                
                // Add boundary source
                if (include_boundary_source_)
                {
                    int index = g + number_of_groups * o;
                    local_sum += source->data()[index];
                }

                // Add local sum into full normal sum
                sum[surface_dimension] += normal * local_sum;
            }
        }
        
        // Add dot product of sum and direction into value
        for (int d = 0; d < dimension; ++d)
        {
            value -= sum[d] * direction[d];
        }
    }
    
    // Add internal source (given contribution)
    int index = g + number_of_groups * (o + number_of_ordinates * i);
    value += x[index];
}

void Weak_RBF_Sweep::
get_matrix_row(int i, // weight function index (row)
               int o, // ordinate
               int g, // group
               vector<int> &indices, // column indices (global basis)
               vector<double> &values) const // column values
{
    // Get data
    shared_ptr<Weight_Function> weight = spatial_discretization_->weight(i);
    vector<double> const is_b_w = weight->is_b_w();
    vector<double> const iv_b_w = weight->iv_b_w();
    vector<double> const iv_b_dw = weight->iv_b_dw();
    vector<double> const iv_db_dw = weight->iv_db_dw();
    vector<double> const direction = angular_discretization_->direction(o);
    vector<double> const sigma_t = weight->material()->sigma_t()->data();
    int number_of_dimensional_moments = spatial_discretization_->number_of_dimensional_moments();
    int number_of_basis_functions = weight->number_of_basis_functions();
    int number_of_boundary_surfaces = weight->number_of_boundary_surfaces();
    int dimension = spatial_discretization_->dimension();
    Weight_Function::Material_Options options = weight->material_options();
    bool include_supg = options.include_supg;
    double tau = options.tau;
    Assert(options.total == Weight_Function::Material_Options::Total::ISOTROPIC); // moment method not yet implemented
    
    // Get indices
    {
        vector<int> const indices_data = weight->basis_function_indices();
        indices = indices_data;
    }
    
    // Get values
    values.assign(number_of_basis_functions, 0);
    for (int j = 0; j < number_of_basis_functions; ++j) // basis function index
    {
        double &value = values[j];
        
        // Add streaming surface contribution
        {
            // Get sum of normals and integrals
            vector<double> sum(dimension, 0);
            for (int s = 0; s < number_of_boundary_surfaces; ++s)
            {
                shared_ptr<Cartesian_Plane> surface = weight->boundary_surface(s);
                int surface_dimension = surface->surface_dimension();
                double const normal = surface->normal();

                // Only for outgoing surfaces
                double dot = normal * direction[surface_dimension];
                if (dot > 0)
                {
                    int is_index = s + number_of_boundary_surfaces * j;
                    sum[surface_dimension] += normal * is_b_w[is_index];
                }
            }
            
            // Add dot product with direction into value
            for (int d = 0; d < dimension; ++d)
            {
                value += sum[d] * direction[d];
            }
        }
        
        // Add streaming volume contribution
        for (int d = 0; d < dimension; ++d)
        {
            int iv_index = d + dimension * j;
            value -= direction[d] * iv_b_dw[iv_index];
        }
        
        // Add streaming SUPG contribution
        if (include_supg)
        {
            for (int d1 = 0; d1 < dimension; ++d1)
            {
                double sum = 0;
                for (int d2 = 0; d2 < dimension; ++d2)
                {
                    int iv_index = d2 + dimension * (d1 + dimension * j);
                    sum += iv_db_dw[iv_index] * direction[d2];
                }
                value += tau * sum * direction[d1];
            }
        }

        // Add absorption term
        {
            int sig_index = 0 + number_of_dimensional_moments * g;
            value += iv_b_w[j] * sigma_t[sig_index];
        }

        // Add absorption SUPG term
        if (include_supg)
        {
            for (int d = 0; d < dimension; ++d)
            {
                int sig_index = (d + 1) + number_of_dimensional_moments * g;
                int iv_index = d + dimension * j;
                value += tau * direction[d] * sigma_t[sig_index] * iv_b_dw[iv_index];
            }
        }
    }
}

void Weak_RBF_Sweep::
output(XML_Node output_node) const
{
    
}

Weak_RBF_Sweep::Sweep_Solver::
Sweep_Solver(Weak_RBF_Sweep const &wrs):
    wrs_(wrs)
{
}

Weak_RBF_Sweep::Trilinos_Solver::
Trilinos_Solver(Weak_RBF_Sweep const &wrs):
    Sweep_Solver(wrs)
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    comm_ = make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);
    map_ = make_shared<Epetra_Map>(number_of_points, 0, *comm_);
    lhs_ = make_shared<Epetra_Vector>(*map_);
    rhs_ = make_shared<Epetra_Vector>(*map_);
    lhs_->PutScalar(1.0);
    rhs_->PutScalar(1.0);
}

shared_ptr<Epetra_CrsMatrix> Weak_RBF_Sweep::Trilinos_Solver::
get_matrix(int o,
           int g) const
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    vector<int> const number_of_basis_functions = wrs_.spatial_discretization_->number_of_basis_functions();
    
    shared_ptr<Epetra_CrsMatrix> mat
        = make_shared<Epetra_CrsMatrix>(Copy, // Data access
                                        *map_,
                                        &number_of_basis_functions[0], // Num entries per row
                                        true); // Static profile
    for (int i = 0; i < number_of_points; ++i)
    {
        vector<int> indices;
        vector<double> values;
        wrs_.get_matrix_row(i,
                            o,
                            g,
                            indices,
                            values);
        mat->InsertGlobalValues(i, // Row
                                number_of_basis_functions[i], // Num entries
                                &values[0],
                                &indices[0]);
    }
    mat->FillComplete();
    mat->OptimizeStorage();

    return mat;
}

void Weak_RBF_Sweep::Trilinos_Solver::
set_rhs(int o,
        int g,
        vector<double> const &x) const
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();

    for (int i = 0; i < number_of_points; ++i)
    {
        double value;
        wrs_.get_rhs(i,
                     o,
                     g,
                     x,
                     value);
        
        (*rhs_)[i] = value;
    }
}

Weak_RBF_Sweep::Amesos_Solver::
Amesos_Solver(Weak_RBF_Sweep const &wrs):
    Trilinos_Solver(wrs)
{
    // Initialize matrices
    int number_of_groups = wrs_.energy_discretization_->number_of_groups();
    int number_of_ordinates = wrs_.angular_discretization_->number_of_ordinates();
    mat_.resize(number_of_groups * number_of_ordinates);
    problem_.resize(number_of_groups * number_of_ordinates);
    solver_.resize(number_of_groups * number_of_ordinates);
    Amesos factory;
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k = g + number_of_groups * o;

            mat_[k] = get_matrix(o,
                                 g);
            problem_[k]
                = make_shared<Epetra_LinearProblem>(mat_[k].get(),
                                                    lhs_.get(),
                                                    rhs_.get());
            solver_[k]
                = make_shared<Amesos_BaseSolver*>(factory.Create("Klu",
                                                                 *problem_[k]));
            (*solver_[k])->SymbolicFactorization();
            (*solver_[k])->NumericFactorization();
        }
    }
}

void Weak_RBF_Sweep::Amesos_Solver::
solve(vector<double> &x) const
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    int number_of_groups = wrs_.energy_discretization_->number_of_groups();
    int number_of_ordinates = wrs_.angular_discretization_->number_of_ordinates();

    // Solve independently for each ordinate and group
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            // Set current RHS value
            set_rhs(o,
                    g,
                    x);
            // Solve, putting result into LHS
            int k = g + number_of_groups * o;
            (*solver_[k])->Solve();
            
            // Update solution value (overwrite x for this o and g)
            for (int i = 0; i < number_of_points; ++i)
            {
                int k_x = g + number_of_groups * (o + number_of_ordinates * i);
                x[k_x] = (*lhs_)[i];
            }
        }
    }
}

Weak_RBF_Sweep::Aztec_Solver::
Aztec_Solver(Weak_RBF_Sweep const &wrs,
             Options options):
    Trilinos_Solver(wrs),
    options_(options)
{
}

void Weak_RBF_Sweep::Aztec_Solver::
solve(vector<double> &x) const
{
    int number_of_points = wrs_.spatial_discretization_->number_of_points();
    int number_of_groups = wrs_.energy_discretization_->number_of_groups();
    int number_of_ordinates = wrs_.angular_discretization_->number_of_ordinates();

    // Solve independently for each ordinate and group
    for (int o = 0; o < number_of_ordinates; ++o)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            // Set current RHS value
            set_rhs(o,
                    g,
                    x);

            // Get matrix
            std::shared_ptr<Epetra_CrsMatrix> mat = get_matrix(o,
                                                               g);

            // Get linear problem
            std::shared_ptr<Epetra_LinearProblem> problem
                = make_shared<Epetra_LinearProblem>(mat.get(),
                                                    lhs_.get(),
                                                    rhs_.get());
            
            // Get solver
            shared_ptr<AztecOO> solver
                = make_shared<AztecOO>(*problem);
            solver->SetAztecOption(AZ_solver, AZ_gmres);
            solver->SetAztecOption(AZ_kspace, options_.kspace);
            solver->SetAztecOption(AZ_precond, AZ_none);
            solver->SetAztecOption(AZ_output, AZ_none);
            
            // Solve, putting result into LHS
            solver->Iterate(options_.max_iterations,
                            options_.tolerance);
            
            // Update solution value (overwrite x for this o and g)
            for (int i = 0; i < number_of_points; ++i)
            {
                int k_x = g + number_of_groups * (o + number_of_ordinates * i);
                x[k_x] = (*lhs_)[i];
            }
        }
    }
}
