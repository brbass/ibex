#include "Weak_RBF_Sweep.hh"

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "Transport_Discretization.hh"
#include "Vector_Funcitons.hh"

namespace vf = Vector_Functions;

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
}

void Weak_RBF_Sweep::
check_class_invariants() const
{
    
}

void Weak_RBF_Sweep::
apply(vector<double> &x) const
{
    
}

void Weak_RBF_Sweep::
get_matrix_row(int i,
               int o,
               int g,
               vector<int> &indices,
               vector<double> &values) const
{
    // Get data
    shared_ptr<Weight_Function> weight = spatial_discretization_->weight(i);
    vector<double> is_b_w = weight->is_b_w();
    vector<double> iv_b_w = weight->iv_b_w();
    vector<double> iv_b_dw = weight->iv_b_dw();
    vector<double> iv_db_dw = weight->iv_db_dw();
    vector<double> direction = angular_discretization_->direction(o);
    vector<double> sigma_t = weight->material()->sigma_t()->data();
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
        vector<double> const indices_data = weight->basis_function_indices();
        indices = indices_data;
    }
    
    // Get values
    values.assign(number_of_basis_functions, 0);
    for (int j = 0; j < number_of_basis_functions; ++j)
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

                int is_index = s + number_of_boundary_surfaces * j;
                sum[surface_dimension] += normal * is_b_w[is_index];
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

