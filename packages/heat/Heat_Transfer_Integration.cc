#include "Heat_Transfer_Integration.hh"

#include <algorithm>
#include <cmath>
#include <limits>

#include "Check.hh"
#include "Heat_Transfer_Data.hh"
#include "Integration_Mesh.hh"
#include "Quadrature_Rule.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weight_Function.hh"

using namespace std;

Heat_Transfer_Integration::
Heat_Transfer_Integration(shared_ptr<Heat_Transfer_Integration_Options> options,
                          shared_ptr<Heat_Transfer_Data> data,
                          shared_ptr<Weak_Spatial_Discretization> spatial):
    options_(options),
    data_(data),
    spatial_(spatial)
{
    // Get integration mesh
    shared_ptr<Integration_Mesh_Options> integration_options
        = make_shared<Integration_Mesh_Options>();
    integration_options->initialize_from_weak_options(spatial->options());
    mesh_ = make_shared<Integration_Mesh>(spatial->dimension(),
                                          spatial->number_of_points(),
                                          integration_options,
                                          spatial->bases(),
                                          spatial->weights());
    Assert(options_);
    Assert(data_);
    Assert(spatial_);
    Assert(mesh_);
    
    initialize_integrals();
    perform_integration();
}

void Heat_Transfer_Integration::
initialize_integrals()
{
    int number_of_points = spatial_->number_of_points();

    matrix_.resize(number_of_points);
    rhs_.assign(number_of_points, 0);
    for (int i = 0; i < number_of_points; ++i)
    {
        shared_ptr<Weight_Function> weight = spatial_->weight(i);
        matrix_[i].assign(weight->number_of_basis_functions(), 0);
    }
}

void Heat_Transfer_Integration::
perform_integration()
{
    int dimension = mesh_->dimension();

    // Perform volume integration
    int number_of_cells = mesh_->number_of_cells();
    for (int i = 0; i < number_of_cells; ++i)
    {
        // Get cell
        shared_ptr<Integration_Cell> const cell = mesh_->cell(i);
        
        // Get quadrature
        int number_of_ordinates;
        vector<vector<double> > ordinates;
        vector<double> weights;
        mesh_->get_volume_quadrature(i,
                                     number_of_ordinates,
                                     ordinates,
                                     weights);
        
        // Get connectivity information
        vector<vector<int> > weight_basis_indices;
        mesh_->get_cell_basis_indices(cell,
                                      weight_basis_indices);

        // Get center positions
        vector<vector<double> > weight_centers;
        vector<vector<double> > basis_centers;
        mesh_->get_basis_weight_centers(cell,
                                        basis_centers,
                                        weight_centers);

        for (int q = 0; q < number_of_ordinates; ++q)
        {
            // Get position
            vector<double> const &position = ordinates[q];
            double const quad_weight = weights[q];

            // For cylindrical, check whether point is inside integration region
            if (options_->geometry == Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_2D)
            {
                double const radius = spatial_->options()->limits[0][1];

                // If point is outside region, continue to next ordinate
                if (position[0] * position[0] + position[1] * position[1] > radius * radius)
                {
                    continue;
                }
            }
            
            // Get values at quadrature point
            vector<double> b_val;
            vector<vector<double> > b_grad;
            vector<double> w_val;
            vector<vector<double> > w_grad;
            mesh_->get_volume_values(cell,
                                     position,
                                     basis_centers,
                                     weight_centers,
                                     b_val,
                                     b_grad,
                                     w_val,
                                     w_grad);
            double const conduction = data_->conduction(position);
            double const source = data_->source(position);

            // Add integrals for each weight function in this cell
            for (int w = 0; w < cell->number_of_weight_functions; ++w)
            {
                // Get global weight function index
                int w_ind = cell->weight_indices[w];

                // Add source to rhs
                switch (options_->geometry)
                {
                case Heat_Transfer_Integration_Options::Geometry::CARTESIAN:
                case Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_2D:
                    rhs_[w_ind] += quad_weight * w_val[w] * source;
                    break;
                case Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_1D:
                    rhs_[w_ind] += quad_weight * w_val[w] * source * position[0];
                    break;
                }
                
                for (int b = 0; b < cell->number_of_basis_functions; ++b)
                {
                    // Get basis index for this weight function
                    int w_b_ind = weight_basis_indices[w][b];

                    // Add convection term to matrix
                    if (w_b_ind != Weight_Function::Errors::DOES_NOT_EXIST)
                    {
                        switch (options_->geometry)
                        {
                        case Heat_Transfer_Integration_Options::Geometry::CARTESIAN:
                        case Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_2D:
                            for (int d = 0; d < dimension; ++d)
                            {
                                matrix_[w_ind][w_b_ind] += quad_weight * w_grad[w][d] * b_grad[b][d] * conduction;
                            }
                            break;
                        case Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_1D:
                            matrix_[w_ind][w_b_ind] += quad_weight * w_grad[w][0] * b_grad[b][0] * conduction * position[0];
                            break;
                        }
                    }
                }
            }
        }
    }

    int number_of_surfaces = mesh_->number_of_surfaces();
    int first_surface;
    switch (options_->geometry)
    {
    case Heat_Transfer_Integration_Options::Geometry::CARTESIAN:
    case Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_2D:
        first_surface = 0;
        break;
    // Skip interior surface for cylindrical 1D geometry
    case Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_1D:
        Assert(number_of_surfaces == 2);
        first_surface = 1;
        break;
    }

    // Split theta into equal parts
    vector<double> limits_t;
    if (options_->geometry == Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_2D)
    {
        limits_t.resize(number_of_surfaces + 1);
        double h = 2 * M_PI / static_cast<double>(number_of_surfaces);
        for (int i = 0; i < number_of_surfaces + 1; ++i)
        {
            limits_t[i] = i * h;
        }
    }
    
    // Perform surface integration
    for (int i = first_surface; i < number_of_surfaces; ++i)
    {
        // Get surface data
        shared_ptr<Integration_Surface> const surface =
            (options_->geometry == Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_2D
             ? get_cylindrical_surface({limits_t[i], limits_t[i+1]},
                                       spatial_->options()->limits[0][1])
             : mesh_->surface(i));
        
        // Get quadrature
        int number_of_ordinates;
        vector<vector<double> > ordinates;
        vector<double> weights;
        switch (options_->geometry)
        {
        case Heat_Transfer_Integration_Options::Geometry::CARTESIAN:
        case Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_1D:
            
            mesh_->get_surface_quadrature(i,
                                          number_of_ordinates,
                                          ordinates,
                                          weights);
            break;
        case Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_2D:
        {
            number_of_ordinates = spatial_->options()->integration_ordinates;
            vector<double> ordinates_t;
            Quadrature_Rule::cartesian_1d(Quadrature_Rule::Quadrature_Type::GAUSS_LEGENDRE,
                                          number_of_ordinates,
                                          limits_t[i],
                                          limits_t[i+1],
                                          ordinates_t,
                                          weights);
            ordinates.resize(number_of_ordinates);
            double const radius = spatial_->options()->limits[0][1];
            for (int q = 0; q < number_of_ordinates; ++q)
            {
                weights[q] *= radius;
                ordinates[q] = {radius * cos(ordinates_t[q]), radius * sin(ordinates_t[q])};
            }
            break;
        }
        }
        
        // Get connectivity information
        vector<vector<int> > weight_basis_indices;
        mesh_->get_surface_basis_indices(surface,
                                         weight_basis_indices);
        
        // Get centers
        vector<vector<double> > weight_centers;
        vector<vector<double> > basis_centers;
        mesh_->get_basis_weight_centers(surface,
                                        basis_centers,
                                        weight_centers);
        
        for (int q = 0; q < number_of_ordinates; ++q)
        {
            // Get position
            double const quad_weight = weights[q];
            vector<double> const &position = ordinates[q];

            // Get basis/weight values at quadrature point
            vector<double> b_val;
            vector<double> w_val;
            mesh_->get_surface_values(surface,
                                      position,
                                      basis_centers,
                                      weight_centers,
                                      b_val,
                                      w_val);
            double const convection = data_->convection(position);
            double const temp_inf = data_->temperature_inf(position);
            
            // Add integrals for each weight function for this surface
            for (int w = 0; w < surface->number_of_weight_functions; ++w)
            {
                // Get global weight function index
                int w_ind = surface->weight_indices[w];
                
                // Add convection term to the source
                switch (options_->geometry)
                {
                case Heat_Transfer_Integration_Options::Geometry::CARTESIAN:
                case Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_2D:
                    rhs_[w_ind] += quad_weight * w_val[w] * convection * temp_inf;
                    break;
                case Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_1D:
                    rhs_[w_ind] += quad_weight * w_val[w] * convection * temp_inf * position[0];
                    break;
                }
                
                for (int b = 0; b < surface->number_of_basis_functions; ++b)
                {
                    // Get basis function index for this weight function
                    int w_b_ind = weight_basis_indices[w][b];
                    
                    // Add convection term to the matrix
                    if (w_b_ind != Weight_Function::Errors::DOES_NOT_EXIST)
                    {
                        switch (options_->geometry)
                        {
                        case Heat_Transfer_Integration_Options::Geometry::CARTESIAN:
                        case Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_2D:
                            matrix_[w_ind][w_b_ind] += quad_weight * w_val[w] * b_val[b] * convection;
                            break;
                        case Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_1D:
                            matrix_[w_ind][w_b_ind] += quad_weight * w_val[w] * b_val[b] * convection * position[0];
                            break;
                        }
                    }
                }
            }
        }
    }
}

std::shared_ptr<Integration_Surface> Heat_Transfer_Integration::
get_cylindrical_surface(vector<double> limit_t,
                        double radius) const
{
    // Requires identical basis functions to work
    Assert(spatial_->options()->identical_basis_functions == Weak_Spatial_Discretization_Options::Identical_Basis_Functions::TRUE);
    
    // Get basis and weight functions that intersect with each point
    int number_of_test_points = 5;
    vector<int> nearest_weights(number_of_test_points);
    double dt = (limit_t[1] - limit_t[0]) / static_cast<double>(number_of_test_points - 1);
    for (int i = 0; i < number_of_test_points; ++i)
    {
        double t = dt * i;
        vector<double> position = {radius * cos(t), radius * sin(t)};
        int nearest_point = spatial_->nearest_point(position);
        nearest_weights[i] = nearest_point;
    }
    sort(nearest_weights.begin(), nearest_weights.end());
    nearest_weights.erase(unique(nearest_weights.begin(), nearest_weights.end()), nearest_weights.end());

    // Get basis and weight functions
    vector<int> indices;
    for (int i : nearest_weights)
    {
        for (int j : spatial_->weight(i)->basis_function_indices())
        {
            indices.push_back(j);
        }
    }
    sort(indices.begin(), indices.end());
    indices.erase(unique(indices.begin(), indices.end()), indices.end());
    int number_of_indices = indices.size();

    // Create surface
    shared_ptr<Integration_Surface> surface = make_shared<Integration_Surface>();
    surface->number_of_basis_functions = number_of_indices;
    surface->number_of_weight_functions = number_of_indices;
    surface->basis_indices = indices;
    surface->weight_indices = indices;

    return surface;
}
