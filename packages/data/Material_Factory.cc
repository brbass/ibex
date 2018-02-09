#include "Material_Factory.hh"

#include "Angular_Discretization.hh"
#include "Cross_Section.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "XML_Node.hh"

using namespace std;

Material_Factory::
Material_Factory(shared_ptr<Angular_Discretization> angular,
                 shared_ptr<Energy_Discretization> energy):
    angular_(angular),
    energy_(energy)
{
}

shared_ptr<Material> Material_Factory::
get_standard_material(int index,
                      vector<double> const &sigma_t_data,
                      vector<double> const &sigma_s_data,
                      vector<double> const &nu_data,
                      vector<double> const &sigma_f_data,
                      vector<double> const &chi_data,
                      vector<double> const &internal_source_data) const
{
    // Get size info
    int number_of_moments = angular_->number_of_moments();
    int number_of_groups = energy_->number_of_groups();
    
    // Get dependencies for standard cross section
    Cross_Section::Dependencies none_group;
    none_group.energy = Cross_Section::Dependencies::Energy::GROUP;
    Cross_Section::Dependencies moment_group;
    moment_group.angular = Cross_Section::Dependencies::Angular::MOMENTS;
    moment_group.energy = Cross_Section::Dependencies::Energy::GROUP;
    Cross_Section::Dependencies scattering_group2;
    scattering_group2.angular = Cross_Section::Dependencies::Angular::SCATTERING_MOMENTS;
    scattering_group2.energy = Cross_Section::Dependencies::Energy::GROUP_TO_GROUP;

    // Make cross sections
    shared_ptr<Cross_Section> sigma_t
        = make_shared<Cross_Section>(none_group,
                                     angular_,
                                     energy_,
                                     sigma_t_data);
    shared_ptr<Cross_Section> sigma_s
        = make_shared<Cross_Section>(scattering_group2,
                                     angular_,
                                     energy_,
                                     sigma_s_data);
    shared_ptr<Cross_Section> nu
        = make_shared<Cross_Section>(none_group,
                                     angular_,
                                     energy_,
                                     nu_data);
    shared_ptr<Cross_Section> sigma_f
        = make_shared<Cross_Section>(none_group,
                                     angular_,
                                     energy_,
                                     sigma_f_data);
    shared_ptr<Cross_Section> chi
        = make_shared<Cross_Section>(none_group,
                                     angular_,
                                     energy_,
                                     chi_data);
    vector<double> internal_source_data_temp = internal_source_data;
    if (internal_source_data.size() == number_of_groups)
    {
        internal_source_data_temp.resize(number_of_groups * number_of_moments, 0);
    }
    shared_ptr<Cross_Section> internal_source
        = make_shared<Cross_Section>(moment_group,
                                     angular_,
                                     energy_,
                                     internal_source_data_temp);

    // Return material
    return make_shared<Material>(index,
                                 angular_,
                                 energy_,
                                 sigma_t,
                                 sigma_s,
                                 nu,
                                 sigma_f,
                                 chi,
                                 internal_source);
}

shared_ptr<Material> Material_Factory::
get_full_fission_material(int index,
                          vector<double> const &sigma_t_data,
                          vector<double> const &sigma_s_data,
                          vector<double> const &sigma_f_data,
                          vector<double> const &internal_source_data) const
{
    // Get size info
    int number_of_moments = angular_->number_of_moments();
    int number_of_groups = energy_->number_of_groups();
    
    // Get dependencies for standard cross section
    Cross_Section::Dependencies none;
    Cross_Section::Dependencies none_group;
    none_group.energy = Cross_Section::Dependencies::Energy::GROUP;
    Cross_Section::Dependencies none_group2;
    none_group2.energy = Cross_Section::Dependencies::Energy::GROUP_TO_GROUP;
    Cross_Section::Dependencies moment_group;
    moment_group.angular = Cross_Section::Dependencies::Angular::MOMENTS;
    moment_group.energy = Cross_Section::Dependencies::Energy::GROUP;
    Cross_Section::Dependencies scattering_group2;
    scattering_group2.angular = Cross_Section::Dependencies::Angular::SCATTERING_MOMENTS;
    scattering_group2.energy = Cross_Section::Dependencies::Energy::GROUP_TO_GROUP;

    // Make cross sections
    shared_ptr<Cross_Section> sigma_t
        = make_shared<Cross_Section>(none_group,
                                     angular_,
                                     energy_,
                                     sigma_t_data);
    shared_ptr<Cross_Section> sigma_s
        = make_shared<Cross_Section>(scattering_group2,
                                     angular_,
                                     energy_,
                                     sigma_s_data);
    shared_ptr<Cross_Section> sigma_f
        = make_shared<Cross_Section>(none_group2,
                                     angular_,
                                     energy_,
                                     sigma_f_data);
    shared_ptr<Cross_Section> nu
        = make_shared<Cross_Section>(none,
                                     angular_,
                                     energy_,
                                     vector<double>(1, 0.));
    shared_ptr<Cross_Section> chi
        = make_shared<Cross_Section>(none,
                                     angular_,
                                     energy_,
                                     vector<double>(1, 0.));
    vector<double> internal_source_data_temp = internal_source_data;
    if (internal_source_data.size() == number_of_groups)
    {
        internal_source_data_temp.resize(number_of_groups * number_of_moments, 0);
    }
    shared_ptr<Cross_Section> internal_source
        = make_shared<Cross_Section>(moment_group,
                                     angular_,
                                     energy_,
                                     internal_source_data_temp);
    
    // Return material
    return make_shared<Material>(index,
                                 angular_,
                                 energy_,
                                 sigma_t,
                                 sigma_s,
                                 nu,
                                 sigma_f,
                                 chi,
                                 internal_source);
}
