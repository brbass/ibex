#include "Scattering.hh"

#include <iostream>
#include <memory>
#include <vector>

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "Spatial_Discretization.hh"

using namespace std;

Scattering::
Scattering(shared_ptr<Spatial_Discretization> spatial_discretization,
           shared_ptr<Angular_Discretization> angular_discretization,
           shared_ptr<Energy_Discretization> energy_discretization,
           Scattering_Type scattering_type):
    Scattering_Operator(spatial_discretization,
                        angular_discretization,
                        energy_discretization,
                        scattering_type)
{
}

void Scattering::
apply_full(vector<double> &x) const
{
    vector<double> y(x);
    
    int number_of_points = spatial_discretization_->number_of_points();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_scattering_moments = angular_discretization_->number_of_scattering_moments();
    vector<int> const scattering_indices = angular_discretization_->scattering_indices();

    for (int i = 0; i < number_of_points; ++i)
    {
        shared_ptr<Material> material = spatial_discretization_->point(i)->material();

        vector<double> const sigma_s = material->sigma_s();
        
        for (int m = 0; m < number_of_moments; ++m)
        {
            int l = scattering_indices[m];

            for (int gt = 0; gt < number_of_groups; ++gt)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    double sum = 0;
                        
                    for (int gf = 0; gf < number_of_groups; ++gf)
                    {
                        int k_phi_from = n + number_of_nodes * (gf + number_of_groups * (m + number_of_moments * i));
                        int k_sigma = gf + number_of_groups * (gt + number_of_groups * l);
                        
                        sum += sigma_s[k_sigma] * y[k_phi_from];
                    }
                    
                    int k_phi_to = n + number_of_nodes * (gt + number_of_groups * (m + number_of_moments * i));
                    
                    x[k_phi_to] = sum;
                }
            }
        }
    }
}

void Scattering::
apply_coherent(vector<double> &x) const
{
    int number_of_points = spatial_discretization_->number_of_points();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_scattering_moments = angular_discretization_->number_of_scattering_moments();
    vector<int> const scattering_indices = angular_discretization_->scattering_indices();
    
    for (int i = 0; i < number_of_points; ++i)
    {
        vector<double> const sigma_s = spatial_discretization_->point(i)->material()->sigma_s();
        
        for (int m = 0; m < number_of_moments; ++m)
        {
            int l = scattering_indices[m];

            for (int g = 0; g < number_of_groups; ++g)
            {
                int k_sigma = g + number_of_groups * (g + number_of_groups * l);
                
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    double sum = 0;
                    
                    int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    
                    x[k_phi] = sigma_s[k_sigma] * x[k_phi];
                }
            }
        }
    }
}

