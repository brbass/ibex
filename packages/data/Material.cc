#include "Material.hh"

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "XML_Node.hh"

using namespace std;

Material::
Material(int index,
         shared_ptr<Angular_Discretization> angular_discretization,
         shared_ptr<Energy_Discretization> energy_discretization,
         vector<double> const &sigma_t,
         vector<double> const &sigma_s,
         vector<double> const &nu,
         vector<double> const &sigma_f,
         vector<double> const &chi,
         vector<double> const &internal_source):
    index_(index),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    sigma_t_(sigma_t),
    sigma_s_(sigma_s),
    nu_(nu),
    sigma_f_(sigma_f),
    chi_(chi),
    internal_source_(internal_source)
{
    check_class_invariants();
}

void Material::
check_class_invariants() const
{
    int number_of_scattering_moments = angular_discretization_->number_of_scattering_moments();
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_groups = energy_discretization_->number_of_groups();
    
    Assert(angular_discretization_);
    Assert(energy_discretization_);
    Assert(sigma_t_.size() == number_of_groups);
    Assert(sigma_s_.size() == number_of_scattering_moments * number_of_groups * number_of_groups);
    Assert(nu_.size() == number_of_groups);
    Assert(sigma_f_.size() == number_of_groups);
    Assert(chi_.size() == number_of_groups);
    Assert(internal_source_.size() == number_of_groups);
}

void Material::
output(XML_Node output_node) const
{
    output_node.set_attribute(index_, "index");
    output_node.set_child_vector(sigma_t_, "sigma_t", "group");
    output_node.set_child_vector(sigma_s_, "sigma_s", "group_from-group_to-moment");
    output_node.set_child_vector(nu_, "nu", "group");
    output_node.set_child_vector(sigma_f_, "sigma_f", "group");
    output_node.set_child_vector(chi_, "chi", "group");
    output_node.set_child_vector(internal_source_, "internal_source", "group");
}
