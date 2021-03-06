#include "SUPG_Internal_Source_Operator.hh"

#if defined(ENABLE_OPENMP)
    #include <omp.h>
#endif

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Cross_Section.hh"
#include "Dimensional_Moments.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "Point.hh"
#include "Spatial_Discretization.hh"

using namespace std;

SUPG_Internal_Source_Operator::
SUPG_Internal_Source_Operator(shared_ptr<Spatial_Discretization> spatial,
                         shared_ptr<Angular_Discretization> angular,
                         shared_ptr<Energy_Discretization> energy):
    spatial_(spatial),
    angular_(angular),
    energy_(energy)
{
    int number_of_points = spatial_->number_of_points();
    int number_of_moments = angular_->number_of_moments();
    int number_of_groups = energy_->number_of_groups();
    int number_of_dimensional_moments = spatial_->dimensional_moments()->number_of_dimensional_moments();
    int number_of_nodes = spatial_->number_of_nodes();

    int phi_size = number_of_points * number_of_moments * number_of_groups * number_of_nodes;
    column_size_ = phi_size;
    row_size_ = phi_size * number_of_dimensional_moments;
    
    check_class_invariants();
}

void SUPG_Internal_Source_Operator::
check_class_invariants() const
{
    Assert(spatial_);
    Assert(angular_);
    Assert(energy_);

    int number_of_points = spatial_->number_of_points();
    for (int i = 0; i < number_of_points; ++i)
    {
        shared_ptr<Material> material = spatial_->point(i)->material();
        Cross_Section::Dependencies dep = material->internal_source()->dependencies();
        // Make sure angular and energy dependencies for each point are correct
        Assert(dep.angular == Cross_Section::Dependencies::Angular::MOMENTS);
        Assert(dep.energy == Cross_Section::Dependencies::Energy::GROUP);
    }
}

void SUPG_Internal_Source_Operator::
apply(vector<double> &x) const
{
    int number_of_points = spatial_->number_of_points();
    int number_of_nodes = spatial_->number_of_nodes();
    int number_of_moments = angular_->number_of_moments();
    int number_of_groups = energy_->number_of_groups();
    int number_of_dimensional_moments = spatial_->dimensional_moments()->number_of_dimensional_moments();
    
    x.assign(number_of_points * number_of_nodes * number_of_moments * number_of_groups * number_of_dimensional_moments, 0);
    #pragma omp parallel for schedule(dynamic, 10)
    for (int i = 0; i < number_of_points; ++i)
    {
        shared_ptr<Cross_Section> internal_source
            = spatial_->point(i)->material()->internal_source();
        vector<double> const data = internal_source->data();
        Cross_Section::Dependencies deps = internal_source->dependencies();
        
        switch (deps.angular)
        {
        case Cross_Section::Dependencies::Angular::MOMENTS:
            for (int g = 0; g < number_of_groups; ++g)
            {
                for (int d = 0; d < number_of_dimensional_moments; ++d)
                {
                    for (int n = 0; n < number_of_nodes; ++n)
                    {
                        for (int m = 0; m < number_of_moments; ++m)
                        {
                            int k_x = n + number_of_nodes * (d + number_of_dimensional_moments * (g + number_of_groups * (m + number_of_moments * i)));
                            int k_dat = d + number_of_dimensional_moments * (g + number_of_groups * m);
                            x[k_x] = data[k_dat];
                        }
                    }
                }
            }
            break;
        default:
            AssertMsg(false, "internal source with given angular dependency not implemented");
            break;
        }
    }
}
