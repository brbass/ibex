#include "Integral_Error_Operator.hh"

#if defined(ENABLE_OPENMP)
    #include <omp.h>
#endif

#include <cmath>
#include <iostream>

#include "Angular_Discretization.hh"
#include "Basis_Function.hh"
#include "Conversion.hh"
#include "Energy_Discretization.hh"
#include "Integration_Mesh.hh"
#include "Weak_Spatial_Discretization.hh"

using std::make_shared;
using std::pair;
using std::shared_ptr;
using std::string;
using std::vector;

Integral_Error_Operator::
Integral_Error_Operator(Options options,
                        shared_ptr<Integration_Mesh_Options> integration_options,
                        shared_ptr<Weak_Spatial_Discretization> spatial,
                        shared_ptr<Angular_Discretization> angular,
                        shared_ptr<Energy_Discretization> energy,
                        std::function<vector<double>(vector<double> const&)> const& solution):
    initialized_(false),
    options_(options),
    integration_options_(integration_options),
    spatial_(spatial),
    angular_(angular),
    energy_(energy),
    solution_(solution)
{
    switch (options_.angular)
    {
    case Options::Angular::NONE:
        angular_size_ = 1;
    case Options::Angular::MOMENTS:
        angular_size_ = angular->number_of_moments();
        break;
    case Options::Angular::ORDINATES:
        angular_size_ = angular->number_of_ordinates();
        break;
    }

    switch (options_.energy)
    {
    case Options::Energy::NONE:
        energy_size_ = 1;
        break;
    case Options::Energy::GROUP:
        energy_size_ = energy->number_of_groups();
        break;
    }

    number_per_point_ = angular_size_ * energy_size_;
    
    column_size_ = (spatial->number_of_points()
                    * number_per_point_);
    int number_of_cells = 1;
    for (int d = 0; d < spatial->dimension(); ++d)
    {
        number_of_cells *= integration_options->dimensional_cells[d];
    }
    row_size_ = (number_of_cells
                 * number_per_point_);
    
    check_class_invariants();
}

void Integral_Error_Operator::
check_class_invariants() const
{
    Assert(integration_options_);
    Assert(spatial_);
    Assert(angular_);
    Assert(energy_);
    Assert(solution_);

    if (initialized_)
    {
        Assert(mesh_);
    }
}

void Integral_Error_Operator::
apply(vector<double> &x) const
{
    // Get size data
    int number_of_points = spatial_->number_of_points();
    
    // Initialize if applicable
    if (!initialized_)
    {
        mesh_ = make_shared<Integration_Mesh>(spatial_->dimension(),
                                              number_of_points,
                                              integration_options_,
                                              spatial_->bases(),
                                              spatial_->weights());
        initialized_ = true;
        check_class_invariants();
    }

    // Ensure row size initialization is correct
    int number_of_cells = mesh_->number_of_cells();
    Assert(row_size_ == number_of_cells * number_per_point_);
    
    // Apply operator
    vector<double> result(row_size_, 0.);
    vector<double> norm(row_size_, 0.);
    #pragma omp parallel for schedule(dynamic, 10)
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

        // Get center positions
        vector<vector<double> > basis_centers;
        mesh_->get_basis_centers(cell,
                                 basis_centers);
        
        for (int q = 0; q < number_of_ordinates; ++q)
        {
            // Get quadrature position and weight
            double const weight = weights[q];
            vector<double> const &position = ordinates[q];
            
            // Get volume values
            vector<double> b_val;
            mesh_->get_basis_values(cell,
                                    position,
                                    basis_centers,
                                    b_val);
            
            // Get flux
            vector<double> flux;
            get_flux(cell,
                     b_val,
                     x,
                     flux);

            // Get expected solution
            vector<double> expected
                = solution_(position);
            
            // Add to integral
            switch (options_.norm)
            {
            case Options::Norm::INTEGRAL:
                for (int j = 0; j < number_per_point_; ++j)
                {
                    int k_res = j + number_per_point_ * i;
                    int k_flux = j;
                    
                    result[k_res] += (expected[k_flux] - flux[k_flux]) * weight;
                }
                break;
            case Options::Norm::L1:
                for (int j = 0; j < number_per_point_; ++j)
                {
                    int k_res = j + number_per_point_ * i;
                    int k_flux = j;
                    
                    result[k_res] += std::abs(expected[k_flux] - flux[k_flux]) * weight;
                }
                break;
            case Options::Norm::L2:
                for (int j = 0; j < number_per_point_; ++j)
                {
                    int k_res = j + number_per_point_ * i;
                    int k_flux = j;
                    
                    double val = expected[k_flux] - flux[k_flux];
                    result[k_res] += val * val * weight;
                }
                break;
            case Options::Norm::LINF:
                for (int j = 0; j < number_per_point_; ++j)
                {
                    int k_res = j + number_per_point_ * i;
                    int k_flux = j;

                    double val = std::abs(expected[k_flux] - flux[k_flux]);
                    if (val > result[k_res])
                    {
                        result[k_res] = val;
                    }
                }
                break;
            }
        }
        
        // Get the volume (sum of weights)
        double volume = 0;
        for (double weight : weights)
        { 
            volume += weight;
        }
        
        // Normalize the result
        switch (options_.norm)
        {
        case Options::Norm::INTEGRAL: // fallthrough intentional
        case Options::Norm::L1:
            for (int j = 0; j < number_per_point_; ++j)
            {
                int k = j + number_per_point_ * i;
                result[k] /= volume;
            }
            break;
        case Options::Norm::L2:
            for (int j = 0; j < number_per_point_; ++j)
            {
                int k = j + number_per_point_ * i;
                result[k] = sqrt(result[k]) / volume;
            }
            break;
        case Options::Norm::LINF:
            // do nothing
            break;
        }
    }
    
    // Put result into "x"
    x.swap(result);
}

void Integral_Error_Operator::
get_flux(shared_ptr<Integration_Cell> const cell,
         vector<double> const &b_val,
         vector<double> const &coeff,
         vector<double> &flux) const
{
    // Calculate flux
    flux.assign(number_per_point_, 0.);
    for (int i = 0; i < cell->number_of_basis_functions; ++i)
    {
        int const j = cell->basis_indices[i];
        for (int l = 0; l < number_per_point_; ++l)
        {
            int const k_f = l; // Local flux index
            int const k_c = l + number_per_point_ * j; // Global coefficient index
            
            flux[k_f] += b_val[i] * coeff[k_c];
        }
    }
}

shared_ptr<Conversion<Integral_Error_Operator::Options::Norm, string> > Integral_Error_Operator::Options::
norm_conversion() const
{
    vector<pair<Norm, string> > conversions
        = {{Norm::INTEGRAL, "integral"},
           {Norm::L1, "l1"},
           {Norm::L2, "l2"},
           {Norm::LINF, "linf"}};
    return make_shared<Conversion<Norm, string> >(conversions);
}

shared_ptr<Conversion<Integral_Error_Operator::Options::Angular, string> > Integral_Error_Operator::Options::
angular_conversion() const
{
    vector<pair<Angular, string> > conversions
        = {{Angular::NONE, "none"},
           {Angular::MOMENTS, "moments"},
           {Angular::ORDINATES, "ordinates"}};
    return make_shared<Conversion<Angular, string> >(conversions);
}

shared_ptr<Conversion<Integral_Error_Operator::Options::Energy, string> > Integral_Error_Operator::Options::
energy_conversion() const
{
    vector<pair<Energy, string> > conversions
        = {{Energy::NONE, "none"},
           {Energy::GROUP, "group"}};
    return make_shared<Conversion<Energy, string> >(conversions);
}
