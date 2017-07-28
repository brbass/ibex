#include "SUPG_Fission.hh"

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

SUPG_Fission::
SUPG_Fission(shared_ptr<Spatial_Discretization> spatial_discretization,
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

void SUPG_Fission::
check_class_invariants() const
{
    Assert(spatial_discretization_);
    Assert(angular_discretization_);
    Assert(energy_discretization_);

    int number_of_points = spatial_discretization_->number_of_points();
    
    // For now, assume that all points have the same energy dependence
    Cross_Section::Dependencies::Energy energy_dep = spatial_discretization_->point(0)->material()->sigma_f()->dependencies().energy;
    for (int i = 0; i < number_of_points; ++i)
    {
        shared_ptr<Material> material = spatial_discretization_->point(i)->material();
        vector<Cross_Section::Dependencies> deps
            = {material->nu()->dependencies(),
               material->sigma_f()->dependencies(),
               material->chi()->dependencies()};
        
        for (Cross_Section::Dependencies &dep : deps)
        {
            Assert(dep.angular == Cross_Section::Dependencies::Angular::NONE);
        }
        if (energy_dep == Cross_Section::Dependencies::Energy::GROUP_TO_GROUP)
        {
            Assert(deps[1].energy == Cross_Section::Dependencies::Energy::GROUP_TO_GROUP);
            Assert(deps[1].dimensional == Cross_Section::Dependencies::Dimensional::SUPG);
        }
        else
        {
            for (Cross_Section::Dependencies &dep : deps)
            {
                Assert(dep.energy == Cross_Section::Dependencies::Energy::GROUP);
                Assert(dep.dimensional == Cross_Section::Dependencies::Dimensional::SUPG);
            }
        }
    }
}

void SUPG_Fission::
apply_full(vector<double> &x) const
{
    switch (spatial_discretization_->point(0)->material()->sigma_f()->dependencies().energy)
    {
    case Cross_Section::Dependencies::Energy::GROUP:
        group_full(x);
        break;
    case Cross_Section::Dependencies::Energy::GROUP_TO_GROUP:
        group_to_group_full(x);
        break;
    }
}

void SUPG_Fission::
apply_coherent(vector<double> &x) const
{
    switch (spatial_discretization_->point(0)->material()->sigma_f()->dependencies().energy)
    {
    case Cross_Section::Dependencies::Energy::GROUP:
        group_coherent(x);
        break;
    case Cross_Section::Dependencies::Energy::GROUP_TO_GROUP:
        group_to_group_coherent(x);
        break;
    }
}

void SUPG_Fission::
group_to_group_full(vector<double> &x) const
{
    // Get size information
    int number_of_points = spatial_discretization_->number_of_points();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    
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
        shared_ptr<Cross_Section> sigma_f_cs = spatial_discretization_->point(i)->material()->sigma_f();
        vector<double> const sigma_f = sigma_f_cs->data();
        
        int const m = 0;
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
                            int k_sigma = d2 + number_of_dimensional_moments * (gf + number_of_groups * gt);
                            
                            sum += sigma_f[k_sigma] * y[k_phi_from];
                        }

                        int d = dimensional_indices[d1 + number_of_dimensional_moments * d2];
                        int k_phi_to = n + number_of_nodes * (d + number_of_double_dimensional_moments * (gt + number_of_groups * (m + number_of_moments * i)));
                        
                        x[k_phi_to] += sum;
                    }
                }
            }
        }
    }
}

void SUPG_Fission::
group_to_group_coherent(vector<double> &x) const
{
    AssertMsg(false, "coherent not implemented in SUPG_Fission");
}

void SUPG_Fission::
group_full(vector<double> &x) const
{
    AssertMsg(false, "group_full in SUPG_Fission not tested");
    
    int number_of_points = spatial_discretization_->number_of_points();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_dimensional_moments = spatial_discretization_->dimensional_moments()->number_of_dimensional_moments();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    
    {
        int m = 0;
        int d = 0;
        for (int i = 0; i < number_of_points; ++i)
        {
            shared_ptr<Material> material = spatial_discretization_->point(i)->material();
            vector<double> const nu = material->nu()->data();
            vector<double> const sigma_f = material->sigma_f()->data();
            vector<double> const chi = material->chi()->data();
            
            for (int n = 0; n < number_of_nodes; ++n)
            {
                // Calculate fission source
                
                double fission_source = 0;
                    
                for (int g = 0; g < number_of_groups; ++g)
                {
                    int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    int k_xs = d + number_of_dimensional_moments * g;
                    
                    fission_source += nu[k_xs] * sigma_f[k_xs] * x[k_phi];
                }
                    
                // Assign flux back to the appropriate group and zeroth moment
                for (int g = 0; g < number_of_groups; ++g)
                {
                    int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    int k_xs = d + number_of_dimensional_moments * g;
                    
                    x[k_phi] = chi[k_xs] * fission_source;
                }
            }
        }
    }
    
    // Zero out other moments
    for (int i = 0; i < number_of_points; ++i)
    {
        for (int m = 1; m < number_of_moments; ++m)
        {
            for (int g = 0; g < number_of_groups; ++g)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    
                    x[k_phi] = 0;
                }
            }
        }
    }
}

void SUPG_Fission::
group_coherent(vector<double> &x) const
{
    AssertMsg(false, "coherent not implemented in SUPG_Fission");
}
