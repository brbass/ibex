#include "VERA_Transport_Result.hh"

#include <cmath>

#include "Angular_Discretization.hh"
#include "Cross_Section.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "Quadrature_Rule.hh"
#include "Solid_Geometry.hh"
#include "Weak_Spatial_Discretization.hh"
#include "XML_Node.hh"

using namespace std;

VERA_Transport_Result::
VERA_Transport_Result(int heat_dimension,
                      double fuel_radius,
                      double pincell_power,
                      shared_ptr<Solid_Geometry> solid,
                      shared_ptr<Angular_Discretization> angular,
                      shared_ptr<Energy_Discretization> energy,
                      shared_ptr<Weak_Spatial_Discretization> spatial,
                      shared_ptr<Solver> solver,
                      shared_ptr<Solver::Result> result):
    heat_dimension_(heat_dimension),
    fuel_radius_(fuel_radius),
    pincell_power_(pincell_power),
    solid_(solid),
    angular_(angular),
    energy_(energy),
    spatial_(spatial),
    solver_(solver),
    result_(result)
{
    // fuel_radius_ = 0.4096;
    if (heat_dimension_ == 1)
    {
        number_of_ordinates_ = 4;
        Quadrature_Rule::cartesian_1d(Quadrature_Rule::Quadrature_Type::GAUSS_LEGENDRE,
                                      number_of_ordinates_,
                                      0,
                                      M_PI / 4,
                                      ordinates_,
                                      weights_);
    }
    
    normalize();
}

double VERA_Transport_Result::
get_fission_energy(vector<double> const &position)
{
    Assert(heat_dimension_ == 2);
    
    // Check whether position is inside the fuel
    double radius2 = position[0] * position[0] + position[1] * position[1];
    if (radius2 > fuel_radius_ * fuel_radius_)
    {
        return 0;
    }

    // Get size information
    int number_of_groups = energy_->number_of_groups();
    int number_of_moments = angular_->number_of_moments();

    // Get angle-independent material at this radius
    shared_ptr<Material> const material
        = solid_->material(position);
    vector<double> const chi_nu_sigma_f
        = material->sigma_f()->data();
    vector<double> kappa_sigma_f(number_of_groups, 0);
    for (int gf = 0; gf < number_of_groups; ++gf)
    {
        for (int gt = 0; gt < number_of_groups; ++gt)
        {
            kappa_sigma_f[gf] += chi_nu_sigma_f[gf + number_of_groups * gt];
        }
            
        kappa_sigma_f[gf] *= kappa_[gf] / nu_[gf] * mev_to_joule_;
    }
    
    // Get flux values
    vector<double> const flux
        = spatial_->expansion_values(number_of_groups * number_of_moments,
                                     position,
                                     result_->coefficients);
    
    // Get fission source values
    double source = 0;
    int const m = 0;
    for (int g = 0; g < number_of_groups; ++g)
    {
        int const k_flux = g + number_of_groups * m;
        source += flux[k_flux] * kappa_sigma_f[g];
    }
    
    return source;
}

double VERA_Transport_Result::
get_radial_fission_energy(double radius)
{
    Assert(heat_dimension_ == 1);
    
    // Check whether position is inside the fuel
    if (radius > fuel_radius_)
    {
        return 0;
    }
        
    // Get size information
    int number_of_groups = energy_->number_of_groups();
    int number_of_moments = angular_->number_of_moments();

    // Get angle-independent material at this radius
    shared_ptr<Material> const material
        = solid_->material({radius, 0});
    vector<double> const chi_nu_sigma_f
        = material->sigma_f()->data();
    vector<double> kappa_sigma_f(number_of_groups, 0);
    for (int gf = 0; gf < number_of_groups; ++gf)
    {
        for (int gt = 0; gt < number_of_groups; ++gt)
        {
            kappa_sigma_f[gf] += chi_nu_sigma_f[gf + number_of_groups * gt];
        }
            
        kappa_sigma_f[gf] *= kappa_[gf] / nu_[gf] * mev_to_joule_;
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
    double norm = 0;
    switch (heat_dimension_)
    {
    case 1:
    {
        // Integrate source radially
        int number_of_ordinates = 256;
        vector<double> ordinates;
        vector<double> weights;
        Quadrature_Rule::cartesian_1d(Quadrature_Rule::Quadrature_Type::GAUSS_LEGENDRE,
                                      number_of_ordinates,
                                      0, // start position r
                                      fuel_radius_,
                                      ordinates,
                                      weights);
        for (int q = 0; q < number_of_ordinates; ++q)
        {
            double radius = ordinates[q];
            double weight = weights[q];
            double source = get_radial_fission_energy(radius);

            norm += source * radius * weight;
        }
        norm *= 2 * M_PI;
        break;
    }
    case 2:
    {
        // Integrate cylindrical source
        int number_of_ordinates_1d = 128;
        vector<double> ordinates_x;
        vector<double> ordinates_y;
        vector<double> weights;
        Quadrature_Rule::cylindrical_2d(Quadrature_Rule::Quadrature_Type::GAUSS_LEGENDRE,
                                        Quadrature_Rule::Quadrature_Type::GAUSS_LEGENDRE,
                                        number_of_ordinates_1d,
                                        number_of_ordinates_1d,
                                        0, // center position x
                                        0, // center position y
                                        0, // start position r
                                        fuel_radius_,
                                        0, // start position theta
                                        2 * M_PI, // end position theta
                                        ordinates_x,
                                        ordinates_y,
                                        weights);
        vector<vector<double> > ordinates;
        Quadrature_Rule::convert_to_position_2d(ordinates_x,
                                                ordinates_y,
                                                ordinates);
        
        int number_of_ordinates = ordinates.size();
        for (int q = 0;  q < number_of_ordinates; ++q)
        {
            double source = get_fission_energy(ordinates[q]);
            double weight = weights[q];
            
            norm += source * weight;
        }
        break;
    }
    default:
    {
        AssertMsg(false, "dimension not found");
        break;
    }
    }
    
    // Normalize flux coefficients
    // double const desired_power = 173.884;
    vector<double> &coefficients = result_->coefficients;
    for (double &coefficient : coefficients)
    {
        coefficient *= pincell_power_ / norm;
    }
}

void VERA_Transport_Result::
output_data(XML_Node output_node)
{
    solid_->output(output_node.append_child("solid_geometry"));
    angular_->output(output_node.append_child("angular_discretization"));
    energy_->output(output_node.append_child("energy_discretization"));
    spatial_->output(output_node.append_child("spatial_discretization"));
    solver_->output_result(output_node.append_child("result"),
                           result_);
}
