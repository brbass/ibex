#include "Fission.hh"

#include <memory>
#include <vector>

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Cross_Section.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "Spatial_Discretization.hh"

using namespace std;

Fission::
Fission(shared_ptr<Spatial_Discretization> spatial_discretization,
        shared_ptr<Angular_Discretization> angular_discretization,
        shared_ptr<Energy_Discretization> energy_discretization,
        Scattering_Type scattering_type):
    Scattering_Operator(spatial_discretization,
                        angular_discretization,
                        energy_discretization,
                        scattering_type)
{
    check_class_invariants();
}

void Fission::
check_class_invariants() const
{
    Assert(spatial_discretization_);
    Assert(angular_discretization_);
    Assert(energy_discretization_);
    
    for (int i = 0; i < number_of_points; ++i)
    {
        shared_ptr<Material> material = spatial_discretization_->point(i)->material();
        vector<Cross_Section::Dependencies> deps
            = {material->nu()->dependencies(),
               material->sigma_f()->dependencies(),
               material->chi()->dependencies()};
        
        for (Cross_Section::Dependencies &dep : deps)
        {
            Assert(dep.angular == Angular::NONE);
            Assert(dep.energy == Angular::GROUP);
        }
    }
}

void Fission::
apply_full(vector<double> &x) const
{
    int number_of_points = spatial_discretization_->number_of_points();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_dimensional_moments = spatial_discretization_->number_of_dimensional_moments();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    
    {
        int m = 0;
        for (int i = 0; i < number_of_points; ++i)
        {
            shared_ptr<Material> material = spatial_discretization_->point(i)->material();
            vector<double> const nu = material->nu()->data();
            vector<double> const sigma_f = material->sigma_f()->data();
            vector<double> const chi = material->chi()->data();
            
            for (int n = 0; n < number_of_nodes; ++n)
            {
                for (int d = 0; d < number_of_dimensional_moments; ++d)
                {
                    // Calculate fission source
                
                    double fission_source = 0;
                    
                    for (int g = 0; g < number_of_groups; ++g)
                    {
                        int k_phi = n + number_of_nodes * (d + number_of_dimensional_moments * (g + number_of_groups * (m + number_of_moments * i)));
                        int k_xs = d + number_of_dimensional_moments * g;
                        
                        fission_source += nu[k_xs] * sigma_f[g] * x[k_phi];
                    }
                    
                    // Assign flux back to the appropriate group and zeroth moment
                    for (int g = 0; g < number_of_groups; ++g)
                    {
                        int k_phi = n + number_of_nodes * (d + number_of_dimensional_moments * (g + number_of_groups * (m + number_of_moments * i)));
                        int k_xs = d + number_of_dimensional_moments * g;
                        
                        x[k_phi] = chi[k_xs] * fission_source;
                    }
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
                    for (int d = 0; d < number_of_dimensional_moments; ++d)
                    {
                        int k_phi = n + number_of_nodes * (d + number_of_dimensional_moments * (g + number_of_groups * (m + number_of_moments * i)));
                        
                        x[k_phi] = 0;
                    }
                }
            }
        }
    }
}

void Fission::
apply_coherent(vector<double> &x) const
{
    int number_of_points = spatial_discretization_->number_of_points();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_dimensional_moments = spatial_discretization_->number_of_dimensional_moments();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();

    // Apply within-group fission to zeroth moment
    {
        int m = 0;
        for (int i = 0; i < number_of_points; ++i)
        {
            shared_ptr<Material> material = spatial_discretization_->point(i)->material();
            vector<double> const nu = material->nu()->data();
            vector<double> const sigma_f = material->sigma_f()->data();
            vector<double> const chi = material->chi()->data();
            
            for (int g = 0; g < number_of_groups; ++g)
            {
                for (int d = 0; d < number_of_dimensional_moments; ++d)
                {
                    int k_xs = d + number_of_dimensional_moments * g;
                    double cs = chi[k_xs] * nu[k_xs] * sigma_f[k_xs];
                
                    for (int n = 0; n < number_of_nodes; ++n)
                    {
                        int k_phi = n + number_of_nodes * (d + number_of_dimensional_moments * (g + number_of_groups * (m + number_of_moments * i)));
                        
                        x[k_phi] = cs * x[k_phi];
                    }
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
                    for (int d = 0; d < number_of_dimensional_moments; ++d)
                    {
                        int k_phi = n + number_of_nodes * (d + number_of_dimensional_moments * (g + number_of_groups * (m + number_of_moments * i)));
                        
                        x[k_phi] = 0;
                    }
                }
            }
        }
    }
}

