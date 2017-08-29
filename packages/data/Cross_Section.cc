#include "Cross_Section.hh"

#include "Angular_Discretization.hh"
#include "Check.hh"
#include "Energy_Discretization.hh"
#include "Conversion.hh"
#include "XML_Node.hh"

using std::make_shared;
using std::pair;
using std::shared_ptr;
using std::string;
using std::vector;

Cross_Section::
Cross_Section(Dependencies dependencies,
              shared_ptr<Angular_Discretization> angular_discretization,
              shared_ptr<Energy_Discretization> energy_discretization,
              vector<double> const &data):
    dependencies_(dependencies),
    angular_discretization_(angular_discretization),
    energy_discretization_(energy_discretization),
    data_(data)
{
    check_class_invariants();

    // if (zero_out_dimensional)
    // {
    //     int m_size = angular_size();
    //     int g_size = energy_size();
    //     int d_size = dimensional_size();
    //     for (int m = 0; m < m_size; ++m)
    //     {
    //         for (int g = 0; g < g_size; ++g)
    //         {
    //             for (int d = 1; d < d_size; ++d)
    //             {
    //                 int k = d + d_size * (g + g_size * m);

    //                 data_[k] = 0;
    //             }
    //         }
    //     }
    // }
}

int Cross_Section::
size() const
{
    return angular_size() * energy_size();
}

int Cross_Section::
angular_size() const
{
    switch (dependencies_.angular)
    {
    case Dependencies::Angular::NONE:
        return 1;
    case Dependencies::Angular::SCATTERING_MOMENTS:
        return angular_discretization_->number_of_scattering_moments();
    case Dependencies::Angular::MOMENTS:
        return angular_discretization_->number_of_moments();
    case Dependencies::Angular::ORDINATES:
        return angular_discretization_->number_of_ordinates();
    }
}

int Cross_Section::
energy_size() const
{
    int num_groups = energy_discretization_->number_of_groups();

    switch (dependencies_.energy)
    {
    case Dependencies::Energy::NONE:
        return 1;
    case Dependencies::Energy::GROUP:
        return num_groups;
    case Dependencies::Energy::GROUP_TO_GROUP:
        return num_groups * num_groups;
    }
}

int Cross_Section::
dimensional_size() const
{
    int dimension = angular_discretization_->dimension();
    
    switch (dependencies_.dimensional)
    {
    case Dependencies::Dimensional::NONE:
        return 1;
    case Dependencies::Dimensional::SUPG:
        return dimension + 1;
    }
}

int Cross_Section::
spatial_size() const
{
    switch (dependencies_.spatial)
    {
    case Dependencies::Spatial::WEIGHT:
        return 1;
    case Dependencies::Spatial::BASIS:
        return 1;
    case Dependencies::Spatial::BASIS_WEIGHT:
        return dependencies_.number_of_basis_functions;
    }
}

string Cross_Section::
angular_string() const
{
    switch (dependencies_.angular)
    {
    case Dependencies::Angular::NONE:
        return "";
    case Dependencies::Angular::SCATTERING_MOMENTS:
        return "scattering_moment";
    case Dependencies::Angular::MOMENTS:
        return "moment";
    case Dependencies::Angular::ORDINATES:
        return "ordinate";
    }
}

string Cross_Section::
energy_string() const
{
    switch (dependencies_.energy)
    {
    case Dependencies::Energy::NONE:
        return "";
    case Dependencies::Energy::GROUP:
        return "-group";
    case Dependencies::Energy::GROUP_TO_GROUP:
        return "-group_to-group_from";
    }
}

string Cross_Section::
dimensional_string() const
{
    switch (dependencies_.dimensional)
    {
    case Dependencies::Dimensional::NONE:
        return "";
    case Dependencies::Dimensional::SUPG:
        return "-supg";
    }
}

string Cross_Section::
spatial_string() const
{
    switch (dependencies_.spatial)
    {
    case Dependencies::Spatial::WEIGHT:
        return "";
    case Dependencies::Spatial::BASIS:
        return "-basis";
    case Dependencies::Spatial::BASIS_WEIGHT:
        return "-basis_weight";
    }
}


void Cross_Section::
check_class_invariants() const
{
    Assert(angular_discretization_);
    Assert(energy_discretization_);
    Assert(data_.size() == angular_size() * energy_size() * dimensional_size() * spatial_size());
}
    
void Cross_Section::
output(XML_Node output_node) const
{
    string description = angular_string() + energy_string() + dimensional_string() + spatial_string();
    
    output_node.set_vector(data_, description);
}

std::shared_ptr<Conversion<Cross_Section::Dependencies::Angular, string> > Cross_Section::Dependencies::
angular_conversion() const
{
    vector<pair<Angular, string> > conversions
        = {{Angular::NONE, "none"},
           {Angular::SCATTERING_MOMENTS, "scattering_moments"},
           {Angular::MOMENTS, "moments"},
           {Angular::ORDINATES, "ordinates"}};
    return make_shared<Conversion<Angular, string> >(conversions);
}

std::shared_ptr<Conversion<Cross_Section::Dependencies::Energy, string> > Cross_Section::Dependencies::
energy_conversion() const
{
    vector<pair<Energy, string> > conversions
        = {{Energy::NONE, "none"},
           {Energy::GROUP, "group"},
           {Energy::GROUP_TO_GROUP, "group_to_group"}};
    return make_shared<Conversion<Energy, string> >(conversions);
}

std::shared_ptr<Conversion<Cross_Section::Dependencies::Dimensional, string> > Cross_Section::Dependencies::
dimensional_conversion() const
{
    vector<pair<Dimensional, string> > conversions
        = {{Dimensional::NONE, "none"},
           {Dimensional::SUPG, "supg"}};
    return make_shared<Conversion<Dimensional, string> >(conversions);
}

std::shared_ptr<Conversion<Cross_Section::Dependencies::Spatial, string> > Cross_Section::Dependencies::
spatial_conversion() const
{
    vector<pair<Spatial, string> > conversions
        = {{Spatial::WEIGHT, "weight"},
           {Spatial::BASIS, "basis"},
           {Spatial::BASIS_WEIGHT, "basis_weight"}};
    return make_shared<Conversion<Spatial, string> >(conversions);
}

