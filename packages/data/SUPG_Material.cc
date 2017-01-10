#include "SUPG_Material.hh"

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "XML_Node.hh"

using namespace std;

SUPG_Material::
SUPG_Material(int index,
              shared_ptr<Angular_Discretization> angular_discretization,
              shared_ptr<Energy_Discretization> energy_discretization,
              vector<double> const &sigma_t,
              vector<double> const &sigma_s,
              vector<double> const &nu,
              vector<double> const &sigma_f,
              vector<double> const &chi,
              vector<double> const &internal_source):
    Material(index,
             angular_discretization,
             energy_discretization,
             sigma_t,
             sigma_s,
             nu,
             sigma_f,
             chi,
             internal_source)
{
    // check_class_invariants() already exists in Material
}

void SUPG_Material::
check_class_invariants() const
{
    int number_of_moments = angular_discretization_->number_of_moments();
    int number_of_groups = energy_discretization_->number_of_groups();
    int dimension = angular_discretization_->dimension();
    int dp1 = dimension + 1;
    
    Assert(angular_discretization_);
    Assert(energy_discretization_);
    Assert(sigma_t_.size() == number_of_moments * number_of_groups * dp1);
    Assert(sigma_s_.size() == number_of_moments * number_of_groups * number_of_groups * dp1);
    Assert(nu_.size() == number_of_groups * dp1);
    Assert(sigma_f_.size() == number_of_groups * dp1);
    Assert(chi_.size() == number_of_groups * dp1);
    Assert(internal_source_.size() == number_of_groups * dp1);
}

void SUPG_Material::
output(XML_Node output_node) const
{
    output_node.set_attribute(index_, "index");
    output_node.set_child_vector(sigma_t_, "sigma_t", "group-moment-dimension");
    output_node.set_child_vector(sigma_s_, "sigma_s", "group_from-group_to-moment-dimension");
    output_node.set_child_vector(nu_, "nu", "group-dimension");
    output_node.set_child_vector(sigma_f_, "sigma_f", "group-dimension");
    output_node.set_child_vector(chi_, "chi", "group-dimension");
    output_node.set_child_vector(internal_source_, "internal_source", "group-dimension");
}
