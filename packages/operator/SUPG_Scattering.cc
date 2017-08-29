#include "SUPG_Scattering.hh"

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
#include "Spatial_Discretization.hh"

using namespace std;

SUPG_Scattering::
SUPG_Scattering(shared_ptr<Spatial_Discretization> spatial_discretization,
           shared_ptr<Angular_Discretization> angular_discretization,
           shared_ptr<Energy_Discretization> energy_discretization,
           Options options):
    SUPG_Scattering_Operator(spatial_discretization,
                        angular_discretization,
                        energy_discretization,
                        options)
{
    check_class_invariants();
}

void SUPG_Scattering::
check_class_invariants() const
{
    Assert(spatial_discretization_);
    Assert(angular_discretization_);
    Assert(energy_discretization_);

    int number_of_points = spatial_discretization_->number_of_points();
    for (int i = 0; i < number_of_points; ++i)
    {
        shared_ptr<Material> material = spatial_discretization_->point(i)->material();
        Cross_Section::Dependencies dep = material->sigma_s()->dependencies();
        // Make sure angular and energy dependencies for each point are correct
        Assert(dep.angular == Cross_Section::Dependencies::Angular::SCATTERING_MOMENTS
               || dep.angular == Cross_Section::Dependencies::Angular::MOMENTS);
        Assert(dep.energy == Cross_Section::Dependencies::Energy::GROUP_TO_GROUP);
        Assert(dep.dimensional == Cross_Section::Dependencies::Dimensional::SUPG);
    }
}

void SUPG_Scattering::
apply_full(vector<double> &x) const
{
    // Get size information
    int number_of_points = spatial_discretization_->number_of_points();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_scattering_moments = angular_discretization_->number_of_scattering_moments();
    vector<int> const scattering_indices = angular_discretization_->scattering_indices();

    // Get dimensional moments
    shared_ptr<Dimensional_Moments> dimensional_moments = spatial_discretization_->dimensional_moments();
    int number_of_dimensional_moments = dimensional_moments->number_of_dimensional_moments();
    int number_of_double_dimensional_moments = dimensional_moments->number_of_double_dimensional_moments();
    vector<int> const dimensional_indices = dimensional_moments->dimensional_indices();
    
    // Copy source flux
    vector<double> y(x);
    x.assign(number_of_points * number_of_nodes * number_of_groups * number_of_moments * number_of_double_dimensional_moments, 0);
    for (int i = 0; i < number_of_points; ++i)
    {
        // Get cross section information
        shared_ptr<Cross_Section> sigma_s_cs = spatial_discretization_->point(i)->material()->sigma_s();
        vector<double> const sigma_s = sigma_s_cs->data();
        Cross_Section::Dependencies::Angular angular_dep = sigma_s_cs->dependencies().angular;
        
        switch (angular_dep)
        {
        case Cross_Section::Dependencies::Angular::SCATTERING_MOMENTS:
        {
            // Perform scattering
            for (int m = 0; m < number_of_moments; ++m)
            {
                int l = scattering_indices[m];
                
                for (int gt = 0; gt < number_of_groups; ++gt)
                {
                    for (int n = 0; n < number_of_nodes; ++n)
                    {
                        for (int d1 = 0; d1 < number_of_dimensional_moments; ++d1)
                        {
                            for (int d2 = 0; d2 < number_of_dimensional_moments; ++d2)
                            {
                                double sum = 0;
                                
                                for (int gf = 0; gf < number_of_groups; ++gf)
                                {
                                    int k_phi_from = n + number_of_nodes * (d1 + number_of_dimensional_moments * (gf + number_of_groups * (m + number_of_moments * i)));
                                    int k_sigma = d2 + number_of_dimensional_moments * (gf + number_of_groups * (gt + number_of_groups * l));
                                    
                                    sum += sigma_s[k_sigma] * y[k_phi_from];
                                }

                                int d = dimensional_indices[d1 + number_of_dimensional_moments * d2];
                                int k_phi_to = n + number_of_nodes * (d + number_of_double_dimensional_moments * (gt + number_of_groups * (m + number_of_moments * i)));
                                
                                x[k_phi_to] += sum;
                            }
                        }
                    }
                }
            }
            break;
        } // SCATTERING_MOMENTS
        case Cross_Section::Dependencies::Angular::MOMENTS:
        {
            // Perform scattering
            for (int m = 0; m < number_of_moments; ++m)
            {
                for (int gt = 0; gt < number_of_groups; ++gt)
                {
                    for (int n = 0; n < number_of_nodes; ++n)
                    {
                        for (int d1 = 0; d1 < number_of_dimensional_moments; ++d1)
                        {
                            for (int d2 = 0; d2 < number_of_dimensional_moments; ++d2)
                            {
                                double sum = 0;
                                
                                for (int gf = 0; gf < number_of_groups; ++gf)
                                {
                                    int k_phi_from = n + number_of_nodes * (d1 + number_of_dimensional_moments * (gf + number_of_groups * (m + number_of_moments * i)));
                                    int k_sigma = d2 + number_of_dimensional_moments * (gf + number_of_groups * (gt + number_of_groups * m));
                                    
                                    sum += sigma_s[k_sigma] * y[k_phi_from];
                                }

                                int d = dimensional_indices[d1 + number_of_dimensional_moments * d2];
                                int k_phi_to = n + number_of_nodes * (d + number_of_double_dimensional_moments * (gt + number_of_groups * (m + number_of_moments * i)));
                                
                                x[k_phi_to] += sum;
                            }
                        }
                    }
                }
            }
            break;
        } // MOMENTS
        default:
            Assert(false);
            break;
        } // switch (angular_dep)
    } // points
}

void SUPG_Scattering::
apply_coherent(vector<double> &x) const
{
    AssertMsg(false, "coherent scattering not yet implemented in SUPG_Scattering");
}
