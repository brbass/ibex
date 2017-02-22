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
    // Get dimensional dependence from first point
    dimensional_dependence_ = spatial_discretization->point(0)->material()->sigma_f()->dependencies().dimensional;

    // Set number of dimensional moments
    int dimension = spatial_discretization->dimension();
    switch (dimensional_dependence_)
    {
    case Cross_Section::Dependencies::Dimensional::NONE:
        number_of_dimensional_moments_ = 1;
        break;
    case Cross_Section::Dependencies::Dimensional::NONE:
        number_of_dimensional_moments_ = 1 + dimension;
        break;
    }
    
    check_class_invariants();
}

void Fission::
check_class_invariants() const
{
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
            Assert(dep.dimensional == dimensional_dependences_);
        }
    }
}

void Fission::
apply_full(vector<double> &x) const
{
    int number_of_points = spatial_discretization_->number_of_points();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();
    
    {
        int m = 0;
        for (int i = 0; i < number_of_points; ++i)
        {
            shared_ptr<Material> material = spatial_discretization_->point(i)->material();
            shared_ptr<Cross_Section> nu = material->nu();
            shared_ptr<Cross_Section> sigma_f = material->sigma_f();
            shared_ptr<Cross_Section> chi = material->chi();
            
            for (int n = 0; n < number_of_nodes; ++n)
            {
                // Calculate fission source
                
                double fission_source = 0;
                
                for (int g = 0; g < number_of_groups; ++g)
                {
                    int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    
                    fission_source += nu[g] * sigma_f[g] * x[k_phi];
                }
                
                // Assign flux back to the appropriate group and zeroth moment
                for (int g = 0; g < number_of_groups; ++g)
                {
                    int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    
                    x[k_phi] = chi[g] * fission_source;
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

void Fission::
apply_coherent(vector<double> &x) const
{
    int number_of_points = spatial_discretization_->number_of_points();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_ordinates = angular_discretization_->number_of_ordinates();

    // Apply within-group fission to zeroth moment
    {
        int m = 0;
        for (int i = 0; i < number_of_points; ++i)
        {
            shared_ptr<Material> material = spatial_discretization_->point(i)->material();

            vector<double> const chi = material->chi();
            vector<double> const nu = material->nu();
            vector<double> const sigma_f = material->sigma_f();
            
            for (int g = 0; g < number_of_groups; ++g)
            {
                double cs = chi[g] * nu[g] * sigma_f[g];
                
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    int k_phi = n + number_of_nodes * (g + number_of_groups * (m + number_of_moments * i));
                    
                    x[k_phi] = cs * x[k_phi];
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

