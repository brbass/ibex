#include "Full_Scattering.hh"

#if defined(ENABLE_OPENMP)
    #include <omp.h>
#endif

#include <iostream>
#include <memory>
#include <vector>

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Cross_Section.hh"
#include "Dimensional_Moments.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weight_Function.hh"

using namespace std;

Full_Scattering::
Full_Scattering(shared_ptr<Weak_Spatial_Discretization> spatial_discretization,
                shared_ptr<Angular_Discretization> angular_discretization,
                shared_ptr<Energy_Discretization> energy_discretization,
                Options options):
    Full_Scattering_Operator(spatial_discretization,
                             angular_discretization,
                             energy_discretization,
                             options)
{
    check_class_invariants();
}

void Full_Scattering::
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
        Assert(dep.angular == Cross_Section::Dependencies::Angular::SCATTERING_MOMENTS);               
        Assert(dep.energy == Cross_Section::Dependencies::Energy::GROUP_TO_GROUP);
        Assert(dep.spatial == Cross_Section::Dependencies::Spatial::BASIS_WEIGHT);
    }
}

void Full_Scattering::
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
    
    // Copy source flux
    vector<double> y(x);
    x.assign(number_of_points * number_of_nodes * number_of_groups * number_of_moments * number_of_dimensional_moments, 0);
    #pragma omp parallel for
    for (int i = 0; i < number_of_points; ++i)
    {
        // Get cross section information
        shared_ptr<Weight_Function> weight = spatial_discretization_->weight(i);
        shared_ptr<Cross_Section> const sigma_s_cs = weight->material()->sigma_s();
        vector<double> const sigma_s = sigma_s_cs->data();
        int const number_of_basis_functions = weight->number_of_basis_functions();
        vector<int> const basis_function_indices = weight->basis_function_indices();
        
        // Perform scattering
        for (int m = 0; m < number_of_moments; ++m)
        {
            int l = scattering_indices[m];
                
            for (int gt = 0; gt < number_of_groups; ++gt)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    for (int d = 0; d < number_of_dimensional_moments; ++d)
                    {
                        double sum = 0;
                        
                        for (int j = 0; j < number_of_basis_functions; ++j)
                        {
                            int b = basis_function_indices[j];
                            
                            for (int gf = 0; gf < number_of_groups; ++gf)
                            {
                                int k_phi_from = n + number_of_nodes * (gf + number_of_groups * (m + number_of_moments * b));
                                int k_sigma = d + number_of_dimensional_moments * (gf + number_of_groups * (gt + number_of_groups * (l + number_of_scattering_moments * j)));
                                
                                sum += sigma_s[k_sigma] * y[k_phi_from];
                            } // from group
                        } // basis functions
                        
                        int k_phi_to = n + number_of_nodes * (d + number_of_dimensional_moments * (gt + number_of_groups * (m + number_of_moments * i)));
                        
                        x[k_phi_to] = sum;
                    } // dimensional moments
                } // nodes
            } // to group
        } // moments
    } // weight functions
}

void Full_Scattering::
apply_coherent(vector<double> &x) const
{
    AssertMsg(false, "coherent scattering not yet implemented in Full_Scattering");
}
