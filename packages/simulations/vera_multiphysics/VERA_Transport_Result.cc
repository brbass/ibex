#include "VERA_Transport_Result.hh"

#include <cmath>

#include "Angular_Discretization.hh"
#include "Cross_Section.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "Quadrature_Rule.hh"
#include "Solid_Geometry.hh"
#include "Weak_Spatial_Discretization.hh"

using namespace std;

VERA_Transport_Result::
VERA_Transport_Result(shared_ptr<Solid_Geometry> solid,
                      shared_ptr<Angular_Discretization> angular,
                      shared_ptr<Energy_Discretization> energy,
                      shared_ptr<Weak_Spatial_Discretization> spatial,
                      shared_ptr<Solver::Result> result):
    solid_(solid),
    angular_(angular),
    energy_(energy),
    spatial_(spatial),
    result_(result)
{
    fuel_radius_ = 0.4096;
    number_of_ordinates_ = 4;
    Quadrature_Rule::cartesian_1d(Quadrature_Rule::Quadrature_Type::GAUSS_LEGENDRE,
                                  number_of_ordinates_,
                                  0,
                                  M_PI / 4,
                                  ordinates_,
                                  weights_);
    
    normalize();
}

double VERA_Transport_Result::
get_radial_fission_energy(double radius)
{
    // Check whether position is inside the fuel
    if (radius > fuel_radius_)
    {
        return 0;
    }
        
    // Get size information
    int number_of_groups = energy_->number_of_groups();
    int number_of_moments = angular_->number_of_moments();

    // Get angle-independent material at this radius
    double const mev_to_joule = 1.6021766e-13;
    shared_ptr<Material> const material
        = solid_->material({radius, 0});
    vector<double> const chi_nu_sigma_f
        = material->sigma_f()->data();
    vector<double> const nu = {2.65063, 2.43223};
    vector<double> const kappa = {196.155, 193.083};
    vector<double> kappa_sigma_f(number_of_groups, 0);
    for (int gf = 0; gf < number_of_groups; ++gf)
    {
        for (int gt = 0; gt < number_of_groups; ++gt)
        {
            kappa_sigma_f[gf] += chi_nu_sigma_f[gf + number_of_groups * gt];
        }
            
        kappa_sigma_f[gf] *= kappa[gf] / nu[gf] * mev_to_joule;
    }
        
    // Integrate source azimuthally from 0 to pi/4
    double source = 0;
    for (int q = 0; q < number_of_ordinates_; ++q)
    {
        // Get position
        double const theta = ordinates_[q];
        vector<double> const position
            = {radius * cos(theta),
               radius * sin(theta)};
            
        // Get flux values
        vector<double> const flux
            = spatial_->expansion_values(number_of_groups * number_of_moments,
                                         position,
                                         result_->coefficients);
            
        // Get fission source values
        int const m = 0;
        for (int g = 0; g < number_of_groups; ++g)
        {
            int const k_flux = g + number_of_groups * m;
            source += flux[k_flux] * kappa_sigma_f[g] * weights_[q];
        }
    }
        
    // Normalize integral
    source *= 4 / M_PI;
        
    return source;
}

void VERA_Transport_Result::
normalize()
{
    // Integrate source radially
    int number_of_ordinates = 256;
    vector<double> ordinates;
    vector<double> weights;
    Quadrature_Rule::cartesian_1d(Quadrature_Rule::Quadrature_Type::GAUSS_LEGENDRE,
                                  number_of_ordinates,
                                  0,
                                  fuel_radius_,
                                  ordinates,
                                  weights);
    double sum = 0;
    for (int q = 0; q < number_of_ordinates; ++q)
    {
        double radius = ordinates[q];
        double weight = weights[q];
        double source = get_radial_fission_energy(radius);

        sum += source * radius * weight;
    }
    sum *= 2 * M_PI;
        
    // Normalize flux coefficients
    double const desired_power = 173.884;
    vector<double> &coefficients = result_->coefficients;
    for (double &coefficient : coefficients)
    {
        coefficient *= desired_power / sum;
    }
}

