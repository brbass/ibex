#include "Combined_SUPG_Fission.hh"

#include <iostream>
#include <memory>
#include <vector>

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Cross_Section.hh"
#include "Dimensional_Moments.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "Point.hh"
#include "Weak_Spatial_Discretization.hh"

using namespace std;

Combined_SUPG_Fission::
Combined_SUPG_Fission(shared_ptr<Weak_Spatial_Discretization> spatial_discretization,
                      shared_ptr<Angular_Discretization> angular_discretization,
                      shared_ptr<Energy_Discretization> energy_discretization,
                      Options options):
    Combined_SUPG_Operator(spatial_discretization,
                           angular_discretization,
                           energy_discretization,
                           options)
{
    check_class_invariants();
}

void Combined_SUPG_Fission::
check_class_invariants() const
{
    Assert(spatial_discretization_);
    Assert(angular_discretization_);
    Assert(energy_discretization_);

    int number_of_points = spatial_discretization_->number_of_points();
    for (int i = 0; i < number_of_points; ++i)
    {
        shared_ptr<Material> material = spatial_discretization_->point(i)->material();
        Cross_Section::Dependencies dep = material->sigma_f()->dependencies();
        // Make sure angular and energy dependencies for each point are correct
        Assert(dep.angular == Cross_Section::Dependencies::Angular::NONE);
        Assert(dep.energy == Cross_Section::Dependencies::Energy::GROUP_TO_GROUP);
        Assert(dep.dimensional == Cross_Section::Dependencies::Dimensional::SUPG);
    }
}

void Combined_SUPG_Fission::
apply_full(vector<double> &x) const
{
    // Get size information
    int const number_of_points = spatial_discretization_->number_of_points();
    int const number_of_nodes = spatial_discretization_->number_of_nodes();
    int const number_of_groups = energy_discretization_->number_of_groups();
    int const number_of_moments = angular_discretization_->number_of_moments();
    int const number_of_ordinates = angular_discretization_->number_of_ordinates();
    double const angular_normalization = angular_discretization_->angular_normalization();
    
    // Get dimensional moments
    shared_ptr<Dimensional_Moments> const dimensional_moments = spatial_discretization_->dimensional_moments();
    int const number_of_dimensional_moments = dimensional_moments->number_of_dimensional_moments();
    
    // Copy source flux
    vector<double> y(x);
    x.assign(number_of_points * number_of_groups * number_of_ordinates, 0);
    int const m = 0;
    for (int i = 0; i < number_of_points; ++i)
    {
        // Get weight function information
        shared_ptr<Weight_Function> const weight = spatial_discretization_->weight(i);
        double const tau = weight->options()->tau;
        
        // Get cross section information
        shared_ptr<Material> const material = weight->material();
        shared_ptr<Cross_Section> const sigma_f_cs = material->sigma_f();
        shared_ptr<Cross_Section> const norm_cs = material->norm();
        vector<double> const sigma_f = sigma_f_cs->data();
        vector<double> const norm = norm_cs->data();
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            vector<double> const direction = angular_discretization_->direction(o);
            vector<double> const coeffs
                = dimensional_moments->coefficients(tau,
                                                    direction);
            
            for (int gt = 0; gt < number_of_groups; ++gt)
            {
                // Apply fission operator
                double fission = 0;
                for (int gf = 0; gf < number_of_groups; ++gf)
                {
                    // Get summations of dimensional moments
                    double phi = 0;
                    double num = 0;
                    double den = 0;
                    for (int d = 0; d < number_of_dimensional_moments; ++d)
                    {
                        int const k_sf = d + number_of_dimensional_moments * (gf + number_of_groups * gt);
                        int const k_phi = d + number_of_dimensional_moments * (gf + number_of_groups * (m + number_of_moments * i));
                        phi += y[k_phi] * coeffs[d];
                        num += sigma_f[k_sf] * coeffs[d];
                        den += norm[k_phi] * coeffs[d];
                    }
                    
                    fission += num / den * phi;
                }
                
                // Apply angular normalization and assign result
                int const k_psi = gt + number_of_groups * (o + number_of_ordinates * i);
                x[k_psi] = fission / angular_normalization;
            }
        }
    }
}

void Combined_SUPG_Fission::
apply_coherent(vector<double> &x) const
{
    AssertMsg(false, "coherent scattering not yet implemented in Combined_SUPG_Fission");
}
