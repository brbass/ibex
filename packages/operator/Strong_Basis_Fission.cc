#include "Strong_Basis_Fission.hh"

#include <iostream>
#include <memory>
#include <vector>
#if defined(ENABLE_OPENMP)
    #include <omp.h>
#endif

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Cross_Section.hh"
#include "Dimensional_Moments.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weight_Function.hh"

using namespace std;

Strong_Basis_Fission::
Strong_Basis_Fission(shared_ptr<Weak_Spatial_Discretization> spatial_discretization,
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

void Strong_Basis_Fission::
check_class_invariants() const
{
    Assert(spatial_discretization_);
    Assert(angular_discretization_);
    Assert(energy_discretization_);
    Assert(spatial_discretization_->options()->weighting
           == Weak_Spatial_Discretization_Options::Weighting::BASIS);
    Assert(spatial_discretization_->options()->discretization
           == Weak_Spatial_Discretization_Options::Discretization::STRONG);
    
    int number_of_points = spatial_discretization_->number_of_points();
    for (int i = 0; i < number_of_points; ++i)
    {
        shared_ptr<Material> material = spatial_discretization_->point(i)->material();
        Cross_Section::Dependencies dep = material->sigma_s()->dependencies();
        // Make sure angular and energy dependencies for each point are correct
        Assert(dep.angular == Cross_Section::Dependencies::Angular::SCATTERING_MOMENTS);
        Assert(dep.energy == Cross_Section::Dependencies::Energy::GROUP_TO_GROUP);
        Assert(dep.spatial == Cross_Section::Dependencies::Spatial::BASIS);
    }
}

void Strong_Basis_Fission::
apply_full(vector<double> &x) const
{
    // Get size information
    int dimension = spatial_discretization_->dimension();
    int number_of_points = spatial_discretization_->number_of_points();
    int number_of_nodes = spatial_discretization_->number_of_nodes();
    int number_of_groups = energy_discretization_->number_of_groups();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_scattering_moments = angular_discretization_->number_of_scattering_moments();
    vector<int> const scattering_indices = angular_discretization_->scattering_indices();
    
    // Copy source flux
    vector<double> y(x);
    x.assign(number_of_points * number_of_nodes * number_of_groups * number_of_moments, 0);
    #pragma omp parallel for schedule(dynamic, 10)
    for (int i = 0; i < number_of_points; ++i)
    {
        // Get weight function and connectivity information
        shared_ptr<Weight_Function> weight = spatial_discretization_->weight(i);
        int const number_of_basis_functions = weight->number_of_basis_functions();
        vector<int> const basis_function_indices = weight->basis_function_indices();
        Weight_Function::Values const values = weight->values();
        
        // Perform scattering
        int const m = 0;
        for (int gt = 0; gt < number_of_groups; ++gt)
        {
            for (int n = 0; n < number_of_nodes; ++n)
            {
                double sum = 0;
                        
                for (int j = 0; j < number_of_basis_functions; ++j)
                {
                    // Get summation constants and other basis data
                    int b = basis_function_indices[j];
                    double mult = values.v_b[j];
                    
                    // Get cross section information: stored in weight functions but weighted by basis functions
                    vector<double> const sigma_f
                        = spatial_discretization_->weight(b)->material()->sigma_f()->data();
                            
                    for (int gf = 0; gf < number_of_groups; ++gf)
                    {
                        int k_phi_from = n + number_of_nodes * (gf + number_of_groups * (m + number_of_moments * b));
                        int k_sigma = gf + number_of_groups * gt;
                            
                        sum += sigma_f[k_sigma] * mult * y[k_phi_from];
                    } // from group
                } // basis functions
                        
                int k_phi_to = n + number_of_nodes * (gt + number_of_groups * (m + number_of_moments * i));
                
                x[k_phi_to] = sum;
            } // nodes
        } // to group
    } // weight functions
}

void Strong_Basis_Fission::
apply_coherent(vector<double> &x) const
{
    AssertMsg(false, "coherent scattering not yet implemented in Strong_Basis_Fission");
}
