#include "Strong_RBF_Sweep.hh"

using namespace std;

Strong_RBF_Sweep::
Strong_RBF_Sweep(Options options,
                 shared_ptr<Weak_Spatial_Discretization> spatial_discretization,
                 shared_ptr<Angular_Discretization> angular_discretization,
                 shared_ptr<Energy_Discretization> energy_discretization,
                 shared_ptr<Transport_Discretization> transport_discretization):
    Weak_RBF_Sweep(options,
                   spatial_discretization,
                   angular_discretization,
                   energy_discretization,
                   transport_discretization)
{
}

void Strong_RBF_Sweep::
get_matrix_row(int i, // weight function index (row)
               int o, // ordinate
               int g, // group
               vector<int> &indices, // column indices (global basis)
               vector<double> &values) const // column values
{
    // Get data
    shared_ptr<Weight_Function> weight = spatial_discretization_->weight(i);
    int number_of_boundary_surfaces = weight->number_of_basis_functions();
    bool boundary_point = number_of_boundary_surfaces > 0;
    Weight_Function::Values const values = weight->values();
    vector<double> const &v_b = values.v_b;
    vector<double> const &v_db = values.v_db;
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
        int surface_dimension = surface->surface_dimension();
        double const normal = surface->normal();

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
        value += sigma_t_data[g] * v_b[j];
    }
}

void Strong_RBF_Sweep::
get_rhs(int i,
        int o,
        int g,
        vector<double> const &x,
        double &value) const
{
   // Get data
    shared_ptr<Weight_Function> weight = spatial_discretization_->weight(i);
    int number_of_boundary_surfaces = weight->number_of_basis_functions();
    bool boundary_point = number_of_boundary_surfaces > 0;
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    int number_of_groups = energy_discretization_->number_of_groups();
    int dimension = spatial_discretization_->dimension();
    int psi_size = row_size() - number_of_augments;
    
    value = 0;
    if (boundary_point)
    {
        Assert(number_of_boundary_surfaces < 2);
        shared_ptr<Cartesian_Plane> boundary_surface = weight->boundary_surface(0);
        shared_ptr<Boundary_Source> boundary_source = weight->boundary_source(0);
        int surface_dimension = surface->surface_dimension();
        double const normal = surface->normal();
        
        // Only for incoming surfaces
        double dot = normal * direction[surface_dimension];
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
            break;
        }
        // else consider as internal point below
    }
    
    int k = g + number_of_groups * (o + number_of_ordinates * i);
    value = x[k];
}
