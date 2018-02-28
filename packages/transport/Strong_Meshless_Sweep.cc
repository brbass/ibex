#include "Strong_Meshless_Sweep.hh"

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
#include "Weak_Spatial_Discretization.hh"

using namespace std;

Strong_Meshless_Sweep::
Strong_Meshless_Sweep(Options options,
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

void  Strong_Meshless_Sweep::
check_class_invariants() const
{
    Assert(spatial_discretization_);
    Assert(angular_discretization_);
    Assert(energy_discretization_);
    Assert(solver_);
    Assert(spatial_discretization_->options()->discretization
           == Weak_Spatial_Discretization_Options::Discretization::STRONG);
}

void Strong_Meshless_Sweep::
get_matrix_row(int i, // weight function index (row)
               int o, // ordinate
               int g, // group
               vector<int> &indices, // column indices (global basis)
               vector<double> &values) const // column values
{
    // Get data
    shared_ptr<Weight_Function> const weight = spatial_discretization_->weight(i);
    int const number_of_boundary_surfaces = weight->number_of_boundary_surfaces();
    bool const boundary_point = number_of_boundary_surfaces > 0;
    Weight_Function::Values const weight_values = weight->values();
    vector<double> const &v_b = weight_values.v_b;
    vector<double> const &v_db = weight_values.v_db;
    vector<double> const direction = angular_discretization_->direction(o);
    int const number_of_basis_functions = weight->number_of_basis_functions();
    vector<int> const basis_indices = weight->basis_function_indices();
    int const dimension = spatial_discretization_->dimension();
    int const number_of_groups = energy_discretization_->number_of_groups();
    shared_ptr<Material> const material = weight->material();
    shared_ptr<Cross_Section> const sigma_t_cs = material->sigma_t();
    vector<double> const sigma_t_data = sigma_t_cs->data();
    
    // Get indices
    indices = basis_indices;

    // Get values
    values.assign(number_of_basis_functions, 0);
    if (boundary_point)
    {
        Assert(number_of_boundary_surfaces < 2);
        shared_ptr<Cartesian_Plane> boundary_surface = weight->boundary_surface(0);
        int surface_dimension = boundary_surface->surface_dimension();
        double const normal = boundary_surface->normal();

        // Only for incoming surfaces
        double dot = normal * direction[surface_dimension];
        if (dot < 0)
        {
            // Add boundary value at this point
            for (int j = 0; j < number_of_basis_functions; ++j)
            {
                double &value = values[j];

                value = v_b[j];
            }
            
            return;
        }
        // else consider as internal point below
    }

    // Add standard streaming operator
    for (int j = 0; j < number_of_basis_functions; ++j)
    {
        double &value = values[j];

        // Add streaming term
        for (int d = 0; d < dimension; ++d)
        {
            int k = d + dimension * j;
            value += direction[d] * v_db[k];
        }

        // Add collision term
        switch (sigma_t_cs->dependencies().spatial)
        {
        case Cross_Section::Dependencies::Spatial::BASIS:
        {
            int const b = basis_indices[j];
            shared_ptr<Material> const basis_material = spatial_discretization_->weight(b)->material();
            shared_ptr<Cross_Section> const basis_sigma_t_cs = basis_material->sigma_t();
            vector<double> const basis_sigma_t_data = basis_sigma_t_cs->data();
            value += basis_sigma_t_data[g] * v_b[j];
        }
        case Cross_Section::Dependencies::Spatial::WEIGHT:
            value += sigma_t_data[g] * v_b[j];
            break;
        default:
            AssertMsg(false, "weighting method not compatible");
        }
    }
}

void Strong_Meshless_Sweep::
get_rhs(int i,
        int o,
        int g,
        vector<double> const &x,
        double &value) const
{
   // Get data
    shared_ptr<Weight_Function> const weight = spatial_discretization_->weight(i);
    vector<double> const direction = angular_discretization_->direction(o);
    int const number_of_boundary_surfaces = weight->number_of_boundary_surfaces();
    bool const boundary_point = number_of_boundary_surfaces > 0;
    int const number_of_ordinates = angular_discretization_->number_of_ordinates();
    int const number_of_groups = energy_discretization_->number_of_groups();
    int const dimension = spatial_discretization_->dimension();
    int const psi_size = transport_discretization_->psi_size();
    
    value = 0;
    if (boundary_point)
    {
        Assert(number_of_boundary_surfaces < 2);
        shared_ptr<Cartesian_Plane> const boundary_surface = weight->boundary_surface(0);
        shared_ptr<Boundary_Source> const boundary_source = weight->boundary_source(0);
        int const surface_dimension = boundary_surface->surface_dimension();
        double const normal = boundary_surface->normal();
        
        // Only for incoming surfaces
        double const dot = normal * direction[surface_dimension];
        if (dot < 0)
        {
            // Add external source contribution
            if (include_boundary_source_)
            {
                vector<double> const source_data = boundary_source->data();
                int const index = g + number_of_groups * o;
                value += source_data[index];
            }
            
            // Add reflective contribution
            bool has_reflection = transport_discretization_->has_reflection();
            if (has_reflection)
            {
                // Get reflection data
                int const number_of_basis_functions = weight->number_of_basis_functions();
                Weight_Function::Values const values = weight->values();
                vector<double> const &v_b = values.v_b;
                double const alpha = boundary_source->alpha()[g];
                vector<double> normal_vec(dimension, 0);
                normal_vec[surface_dimension] = normal;
                int const o_ref = angular_discretization_->reflect_ordinate(o,
                                                                            normal_vec);
                for (int j = 0; j < number_of_basis_functions; ++j)
                {
                    int const index = psi_size + g + number_of_groups * (o_ref + number_of_ordinates * j);
                    value += alpha * v_b[j] * x[index];
                }
            }
            return;
        }
        // else consider as internal point below
    }
    
    int const k = g + number_of_groups * (o + number_of_ordinates * i);
    value = x[k];
}
