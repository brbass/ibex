#include "Heat_Transfer_Integration.hh"

#include <limits>
#include <iostream>

#include "Check.hh"
#include "Heat_Transfer_Data.hh"
#include "Integration_Mesh.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weight_Function.hh"

using namespace std;

Heat_Transfer_Integration::
Heat_Transfer_Integration(shared_ptr<Heat_Transfer_Integration_Options> options,
                          shared_ptr<Integration_Mesh_Options> integration_options,
                          shared_ptr<Heat_Transfer_Data> data,
                          shared_ptr<Weak_Spatial_Discretization> spatial):
    options_(options),
    mesh_(make_shared<Integration_Mesh>(spatial->dimension(),
                                        spatial->number_of_points(),
                                        integration_options,
                                        spatial->bases(),
                                        spatial->weights())),
    data_(data),
    spatial_(spatial)
{
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
        shared_ptr<Integration_Mesh::Cell> const cell = mesh_->cell(i);
        
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
                int w_ind = cell->weight_indices[i];

                // Add source to rhs
                switch (options_->geometry)
                {
                case Heat_Transfer_Integration_Options::Geometry::CARTESIAN:
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
    
    // Skip interior surface for cylindrical geometry
    int number_of_surfaces = mesh_->number_of_surfaces();
    int first_surface;
    switch (options_->geometry)
    {
    case Heat_Transfer_Integration_Options::Geometry::CARTESIAN:
        first_surface = 0;
        break;
    case Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_1D:
        Assert(number_of_surfaces == 2);
        first_surface = 1;
        break;
    }

    // Perform surface integration
    for (int i = first_surface; i < number_of_surfaces; ++i)
    {
        // Get surface data
        shared_ptr<Integration_Mesh::Surface> const surface = mesh_->surface(i);

        // Get quadrature
        int number_of_ordinates;
        vector<vector<double> > ordinates;
        vector<double> weights;
        mesh_->get_surface_quadrature(i,
                                      number_of_ordinates,
                                      ordinates,
                                      weights);
        
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
                            matrix_[w_ind][w_b_ind] += quad_weight * w_val[w] * b_val[b] * convection;
                            break;
                        case Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_1D:
                            matrix_[w_ind][w_b_ind] += quad_weight * w_val[w] * b_val[b] * position[0];
                            break;
                        }
                    }
                }
            }
        }
    }
}

