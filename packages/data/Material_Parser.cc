#include "Material_Parser.hh"

#include "Angular_Discretization.hh"
#include "Cross_Section.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "XML_Node.hh"

using namespace std;

Material_Parser::
Material_Parser(shared_ptr<Angular_Discretization> angular,
                shared_ptr<Energy_Discretization> energy):
    angular_(angular),
    energy_(energy)
{
}

vector<shared_ptr<Material> > Material_Parser::
parse_from_xml(XML_Node input_node)
{
    // Get size data
    int number_of_materials = input_node.get_child_value<int>("number_of_materials");
    
    int number_of_scattering_moments = angular_->number_of_scattering_moments();
    int number_of_moments = angular_->number_of_moments();
    int number_of_groups = energy_->number_of_groups();

    // Get possible dependencies
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
    
    // Parse the data for each material
    
    vector<shared_ptr<Material> > materials(number_of_materials);
    
    int checksum = 0;
    for (XML_Node material_node = input_node.get_child("material");
         material_node;
         material_node = material_node.get_sibling("material",
                                                   false))
    {
        int a = material_node.get_attribute<int>("index");
        Assert(a < number_of_materials);
        
        // Get cross sections
        shared_ptr<Cross_Section> sigma_t
            = make_shared<Cross_Section>(none_group,
                                         angular_,
                                         energy_,
                                         material_node.get_child_vector<double>("sigma_t", number_of_groups));
        shared_ptr<Cross_Section> sigma_s
            = make_shared<Cross_Section>(scattering_group2,
                                         angular_,
                                         energy_,
                                         material_node.get_child_vector<double>("sigma_s", number_of_groups * number_of_groups * number_of_scattering_moments));
        vector<double> internal_source_data = material_node.get_child_vector<double>("internal_source", number_of_groups);
        internal_source_data.resize(number_of_groups * number_of_moments, 0);
        shared_ptr<Cross_Section> internal_source
            = make_shared<Cross_Section>(moment_group,
                                         angular_,
                                         energy_,
                                         internal_source_data);
        
        // Get fission cross section, checking for full fission matrix
        shared_ptr<Cross_Section> nu;
        shared_ptr<Cross_Section> sigma_f;
        shared_ptr<Cross_Section> chi;
        
        if (material_node.get_child("chi_nu_sigma_f", false))
        {
            nu = make_shared<Cross_Section>(none,
                                            angular_,
                                            energy_,
                                            vector<double>(1, 0.));
            sigma_f = make_shared<Cross_Section>(none_group2,
                                                 angular_,
                                                 energy_,
                                                 material_node.get_child_vector<double>("chi_nu_sigma_f", number_of_groups * number_of_groups));
            chi = make_shared<Cross_Section>(none,
                                             angular_,
                                             energy_,
                                             vector<double>(1, 0.));
        }
        else
        {
            nu = make_shared<Cross_Section>(none_group,
                                            angular_,
                                            energy_,
                                            material_node.get_child_vector<double>("nu", number_of_groups));
            sigma_f = make_shared<Cross_Section>(none_group,
                                                 angular_,
                                                 energy_,
                                                 material_node.get_child_vector<double>("sigma_f", number_of_groups));
            chi = make_shared<Cross_Section>(none_group,
                                             angular_,
                                             energy_,
                                             material_node.get_child_vector<double>("chi", number_of_groups));
        }

        // Create material
        shared_ptr<Material> material = make_shared<Material>(a,
                                                              angular_,
                                                              energy_,
                                                              sigma_t,
                                                              sigma_s,
                                                              nu,
                                                              sigma_f,
                                                              chi,
                                                              internal_source);
        
        materials[a] = material;
        
        checksum += a;
    } // materials
    
    int checksum_expected = number_of_materials * (number_of_materials - 1) / 2;
    AssertMsg(checksum == checksum_expected, "Material indexing incorrect");

    // Check that energy dependencies of sigma_f are the same
    Cross_Section::Dependencies::Energy energy_dep = materials[0]->sigma_f()->dependencies().energy;
    for (shared_ptr<Material> material : materials)
    {
        Assert(material->sigma_f()->dependencies().energy == energy_dep);
    }
    
    return materials;
}
