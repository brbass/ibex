#include "Material.hh"

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Cross_Section.hh"
#include "Energy_Discretization.hh"
#include "XML_Node.hh"

using namespace std;

Material::
Material(int index,
         shared_ptr<Angular_Discretization> angular_discretization,
         shared_ptr<Energy_Discretization> energy_discretization,
         shared_ptr<Cross_Section> sigma_t,
         shared_ptr<Cross_Section> sigma_s,
         shared_ptr<Cross_Section> nu,
         shared_ptr<Cross_Section> sigma_f,
         shared_ptr<Cross_Section> chi,
         shared_ptr<Cross_Section> internal_source,
         shared_ptr<Cross_Section> norm):
    index_(index),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    sigma_t_(sigma_t),
    sigma_s_(sigma_s),
    nu_(nu),
    sigma_f_(sigma_f),
    chi_(chi),
    internal_source_(internal_source),
    norm_(norm)
{
    check_class_invariants();
}

void Material::
check_class_invariants() const
{
    Assert(angular_discretization_);
    Assert(energy_discretization_);
    Assert(sigma_t_);
    Assert(sigma_s_);
    Assert(nu_);
    Assert(sigma_f_);
    Assert(chi_);
    Assert(internal_source_);
}

void Material::
output(XML_Node output_node) const
{
    output_node.set_attribute(index_, "index");
    sigma_t_->output(output_node.append_child("sigma_t"));
    sigma_s_->output(output_node.append_child("sigma_s"));
    nu_->output(output_node.append_child("nu"));
    sigma_f_->output(output_node.append_child("sigma_f"));
    chi_->output(output_node.append_child("chi"));
    internal_source_->output(output_node.append_child("internal_source"));
    if (norm_)
    {
        norm_->output(output_node.append_child("norm"));
    }
}
