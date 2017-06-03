#include "Internal_Source_Operator.hh"

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Cross_Section.hh"
#include "Dimensional_Moments.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "Point.hh"
#include "Spatial_Discretization.hh"

using namespace std;

Internal_Source_Operator::
Internal_Source_Operator(shared_ptr<Spatial_Discretization> spatial,
                         shared_ptr<Angular_Discretization> angular,
                         shared_ptr<Energy_Discretization> energy):
    spatial_(spatial),
    angular_(angular),
    energy_(energy)
{
    int number_of_points = spatial_->number_of_points();
    int number_of_moments = angular_->number_of_moments();
    int number_of_groups = energy_->number_of_groups();
    int number_of_nodes = spatial_->number_of_nodes();
    
    size_ = number_of_points * number_of_moments * number_of_groups * number_of_nodes;

    check_class_invariants();
}

void Internal_Source_Operator::
check_class_invariants() const
{
    Assert(spatial_);
    Assert(angular_);
    Assert(energy_);
}

void Internal_Source_Operator::
apply(vector<double> &x) const
{
    int number_of_points = spatial_->number_of_points();
    int number_of_nodes = spatial_->number_of_nodes();
    int number_of_moments = angular_->number_of_moments();
    int number_of_groups = energy_->number_of_groups();
    int number_of_dimensional_moments = spatial_->dimensional_moments()->number_of_dimensional_moments();
    double normalization = angular_->angular_normalization();
    
    for (int i = 0; i < number_of_points; ++i)
    {
        shared_ptr<Cross_Section> internal_source
            = spatial_->point(i)->material()->internal_source();
        vector<double> const data = internal_source->data();
        Cross_Section::Dependencies deps = internal_source->dependencies();

        int d = 0;
        switch (deps.angular)
        {
        case Cross_Section::Dependencies::Angular::NONE:
            for (int g = 0; g < number_of_groups; ++g)
            {
                for (int n = 0; n < number_of_nodes; ++n)
                {
                    {
                        int m = 0;
                        int k_x = n + number_of_nodes * (d + number_of_dimensional_moments * (g + number_of_groups * (m + number_of_moments * i)));
                        int k_dat = d + number_of_dimensional_moments * g;
                        x[k_x] = data[k_dat];
                    }
                        
                    for (int m = 1; m < number_of_moments; ++m)
                    {
                        int k_x = n + number_of_nodes * (d + number_of_dimensional_moments * (g + number_of_groups * (m + number_of_moments * i)));
                        x[k_x] = 0;
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
