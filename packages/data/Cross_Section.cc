#include "Cross_Section.hh"

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "XML_Node.hh"

using std::shared_ptr;
using std::string;
using std::vector;

Cross_Section::
Cross_Section(Angular_Dependence angular_dependence,
              Energy_Dependence energy_dependence,
              shared_ptr<Angular_Discretization> angular_discretization,
              shared_ptr<Energy_Discretization> energy_discretization,
              vector<double> const &data):
    angular_dependence_(angular_dependence),
    energy_dependence_(energy_dependence),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    data_(data)
{
    check_class_invariants();
}

int Cross_Section::
size() const
{
    return angular_size() * energy_size();
}

int Cross_Section::
angular_size() const
{
    switch (angular_dependence_)
    {
    case Angular_Dependence::NONE:
        return 1;
    case Angular_Dependence::SCATTERING_MOMENTS:
        return angular_discretization_->number_of_scattering_moments();
    case Angular_Dependence::MOMENTS:
        return angular_discretization_->number_of_moments();
    case Angular_Dependence::ORDINATES:
        return angular_discretization_->number_of_ordinates();
    }
}

int Cross_Section::
energy_size() const
{
    int num_groups = energy_discretization_->number_of_groups();

    switch (energy_dependence_)
    {
    case Energy_Dependence::NONE:
        return 1;
    case Energy_Dependence::GROUP:
        return num_groups;
    case Energy_Dependence::GROUP_TO_GROUP:
        return num_groups * num_groups;
    }
}

string Cross_Section::
angular_string() const
{
    switch (angular_dependence_)
    {
    case Angular_Dependence::NONE:
        return "none";
    case Angular_Dependence::SCATTERING_MOMENTS:
        return "scattering_moment";
    case Angular_Dependence::MOMENTS:
        return "moment";
    case Angular_Dependence::ORDINATES:
        return "ordinate";
    }
}

string Cross_Section::
energy_string() const
{
    switch (energy_dependence_)
    {
    case Energy_Dependence::NONE:
        return "none";
    case Energy_Dependence::GROUP:
        return "group";
    case Energy_Dependence::GROUP_TO_GROUP:
        return "group_to-group_from";
    }
}

void Cross_Section::
check_class_invariants() const
{
    Assert(angular_discretization_);
    Assert(energy_discretization_);
    Assert(data_.size() == angular_size() * energy_size());
}
    
void Cross_Section::
output(XML_Node output_node) const
{
    string description = angular_string() + "-" + energy_string();
    
    output_node.set_child_vector(data_, "data", description);
}
