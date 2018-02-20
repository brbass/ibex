#include "Weak_Meshless_Sweep.hh"

#include <iostream>
#if defined(ENABLE_OPENMP)
    #include <omp.h>
#else
    inline int omp_get_num_threads() {return 1;}
    inline int omp_get_thread_num() {return 0;}
#endif

#include "Amesos.h"
#include "AztecOO.h"
#include "AztecOO_ConditionNumber.h"
#include "BelosSolverFactory.hpp"
#include "BelosEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Vector.h"
#include "Ifpack.h"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_RCPStdSharedPtrConversions.hpp"
#include "Teuchos_ArrayRCP.hpp"

#include "Angular_Discretization.hh"
#include "Basis_Function.hh"
#include "Boundary_Source.hh"
#include "Cartesian_Plane.hh"
#include "Check.hh"
#include "Conversion.hh"
#include "Cross_Section.hh"
#include "Dimensional_Moments.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "Transport_Discretization.hh"
#include "XML_Node.hh"

using std::make_shared;
using std::pair;
using std::shared_ptr;
using std::string;
using std::vector;

Weak_Meshless_Sweep::
Weak_Meshless_Sweep(Options options,
               shared_ptr<Weak_Spatial_Discretization> spatial_discretization,
               shared_ptr<Angular_Discretization> angular_discretization,
               shared_ptr<Energy_Discretization> energy_discretization,
               shared_ptr<Transport_Discretization> transport_discretization):
    Meshless_Sweep(options,
                   spatial_discretization,
                   angular_discretization,
                   energy_discretization,
                   transport_discretization)
{
    initialize_solver();
    check_class_invariants();
}

void Weak_Meshless_Sweep::
check_class_invariants() const
{
    Assert(spatial_discretization_);
    Assert(angular_discretization_);
    Assert(energy_discretization_);
    Assert(solver_);
    Assert(spatial_discretization_->options()->discretization
           == Weak_Spatial_Discretization_Options::Discretization::WEAK);
}

void Weak_Meshless_Sweep::
get_rhs(int i,
        int o,
        int g,
        vector<double> const &x,
        double &value) const
{
    // Get data
    shared_ptr<Weight_Function> weight = spatial_discretization_->weight(i);
    vector<double> const &is_b_w = weight->integrals().is_b_w;
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
                shared_ptr<Boundary_Source> source = weight->boundary_source(s);
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

void Weak_Meshless_Sweep::
get_matrix_row(int i, // weight function index (row)
               int o, // ordinate
               int g, // group
               vector<int> &indices, // column indices (global basis)
               vector<double> &values) const // column values
{
    // Get data
    shared_ptr<Weight_Function> weight = spatial_discretization_->weight(i);
    Weight_Function::Integrals const integrals = weight->integrals();
    vector<double> const &iv_w = integrals.iv_w;
    vector<double> const &iv_dw = integrals.iv_dw;
    vector<double> const &is_b_w = integrals.is_b_w;
    vector<double> const &iv_b_w = integrals.iv_b_w;
    vector<double> const &iv_b_dw = integrals.iv_b_dw;
    vector<double> const &iv_db_dw = integrals.iv_db_dw;
    vector<double> const direction = angular_discretization_->direction(o);
    shared_ptr<Dimensional_Moments> const dimensional_moments
        = spatial_discretization_->dimensional_moments();
    int const number_of_dimensional_moments = dimensional_moments->number_of_dimensional_moments();
    int const number_of_basis_functions = weight->number_of_basis_functions();
    vector<int> const basis_indices = weight->basis_function_indices();
    int const number_of_boundary_surfaces = weight->number_of_boundary_surfaces();
    int const dimension = spatial_discretization_->dimension();
    int const number_of_groups = energy_discretization_->number_of_groups();
    shared_ptr<Weak_Spatial_Discretization_Options> const weak_options
        = spatial_discretization_->options();
    shared_ptr<Weight_Function_Options> const weight_options
        = weight->options();
    shared_ptr<Material> const material = weight->material();
    shared_ptr<Cross_Section> const sigma_t_cs = material->sigma_t();
    shared_ptr<Cross_Section> const norm_cs = material->norm();
    vector<double> const sigma_t_data = sigma_t_cs->data();
    
    bool const include_supg = weak_options->include_supg;
    bool const normalized = weak_options->normalized;
    double const tau = weight_options->tau;
    vector<double> const dimensional_coefficients
        = dimensional_moments->coefficients(tau,
                                            direction);

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

        // Add collision term
        switch (sigma_t_cs->dependencies().spatial)
        {
        case Cross_Section::Dependencies::Spatial::BASIS_WEIGHT:
        {
            double sum = 0;
            
            for (int d = 0; d < number_of_dimensional_moments; ++d)
            {
                int k_sigma = d + number_of_dimensional_moments * (g + number_of_groups * j);
                sum += dimensional_coefficients[d] * sigma_t_data[k_sigma];
            }
            
            value += sum;
            break;
        } // Spatial::BASIS_WEIGHT
        case Cross_Section::Dependencies::Spatial::BASIS:
        {
            int const b = basis_indices[j];
            shared_ptr<Material> const basis_material = spatial_discretization_->weight(b)->material();
            shared_ptr<Cross_Section> const basis_sigma_t_cs = basis_material->sigma_t();
            vector<double> const basis_sigma_t_data = basis_sigma_t_cs->data();

            double sum = 0;
            for (int d = 0; d < number_of_dimensional_moments; ++d)
            {
                int const k_sigma = d + number_of_dimensional_moments * g;
                double const mult = (d == 0
                                     ? iv_b_w[j]
                                     : iv_b_dw[d - 1 + dimension * j]);
                sum += dimensional_coefficients[d] * mult * basis_sigma_t_data[k_sigma];
            }
            
            value += sum;
            break;
        }
        case Cross_Section::Dependencies::Spatial::WEIGHT:
        {
            Assert(weak_options->total == Weak_Spatial_Discretization_Options::Total::ISOTROPIC); // moment method not yet implemented
            
            // Get total cross section: leave out higher moments for now
            double sigma_t = 0;
            for (int d = 0; d < number_of_dimensional_moments; ++d)
            {
                int const k_sigma = d + number_of_dimensional_moments * g;
                sigma_t += sigma_t_data[k_sigma] * dimensional_coefficients[d];
            }
            
            // Normalize total cross section if needed
            if (!normalized)
            {
                vector<double> const norm_data = norm_cs->data();
                double norm = 0;
                switch (norm_cs->dependencies().energy)
                {
                case Cross_Section::Dependencies::Energy::NONE:
                    // Norm depends only on dimensional moment
                    for (int d = 0; d < number_of_dimensional_moments; ++d)
                    {
                        norm += norm_data[d] * dimensional_coefficients[d];
                    }
                    break;
                case Cross_Section::Dependencies::Energy::GROUP:
                    // Norm depends on dimensional moment, angular moment and group
                    // Ignore the angular moment for now
                    for (int d = 0; d < number_of_dimensional_moments; ++d)
                    {
                        int const k_norm = d + number_of_dimensional_moments * g;
                        norm += norm_data[k_norm] * dimensional_coefficients[d];
                    }
                    break;
                default:
                    AssertMsg(false, "norm dependency incorrect");
                    break;
                }
                sigma_t /= norm;
            } // if !normalized
    
            double sum = iv_b_w[j];
            if (include_supg)
            {
                for (int d = 0; d < dimension; ++d)
                {
                    int iv_index = d + dimension * j;
                    sum += tau * direction[d] * iv_b_dw[iv_index];
                }
            }
            
            value += sum * sigma_t;
            break;
        } // Spatial::WEIGHT
        } // switch Spatial
    } // basis functions
}

