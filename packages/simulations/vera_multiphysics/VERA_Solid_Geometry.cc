#include "VERA_Solid_Geometry.hh"

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Boundary_Source.hh"
#include "Material.hh"
#include "Material_Factory.hh"

using namespace std;

VERA_Solid_Geometry::
VERA_Solid_Geometry(bool include_ifba,
                    function<double(vector<double> const &)> temperature,
                    shared_ptr<Angular_Discretization> angular,
                    shared_ptr<Energy_Discretization> energy):
    include_ifba_(include_ifba),
    temperature_(temperature),
    angular_(angular),
    energy_(energy),
    material_factory_(make_shared<Material_Factory>(angular,
                                                    energy))
{
    initialize_materials();

    check_class_invariants();
}

int VERA_Solid_Geometry::
dimension() const
{
    return angular_->dimension();
}

void VERA_Solid_Geometry::
check_class_invariants() const
{
    Assert(angular_->dimension() == 2);
    Assert(angular_->number_of_scattering_moments() == 2);
    Assert(energy_->number_of_groups() == 2);
}

void VERA_Solid_Geometry::
initialize_materials()
{
    number_of_material_types_ = 5;
    number_of_xs_types_ = 4;
    materials_.resize(number_of_material_types_,
                      vector<shared_ptr<Material> >(number_of_xs_types_));
    materials_[FUEL][NORM_600]
        = material_factory_->get_full_fission_material(
            0, // index
            {0.3996976, 0.581826884},
            {0.383829618, 0.0,
                    0.000830382, 0.405420306,
                    0.049482651, 0.0,
                    -0.000261476, 0.006013128},
            {0.013687612, 0.255838918,
                    0.000000015, 0.000000189},
            {0, 0});
    
            
                    
}

void VERA_Solid_Geometry::
get_cross_sections(double temperature,
                   vector<double> const &position,
                   vector<double> &sigma_t,
                   vector<double> &sigma_s,
                   vector<double> &sigma_f) const
{
    
}

shared_ptr<Material> VERA_Solid_Geometry::
material(vector<double> const &position) const
{
    int number_of_scattering_moments = angular_->number_of_scattering_moments();
    int number_of_groups = energy_->number_of_groups();
    
    // Initialize data arrays
    vector<double> sigma_t_d(number_of_groups);
    vector<double> sigma_s_d(number_of_groups * number_of_groups * number_of_scattering_moments);
    vector<double> sigma_f_d(number_of_groups * number_of_groups);
    
    // Get cross sections
    double temperature = temperature_(position);
    get_cross_sections(temperature,
                       position,
                       sigma_t_d,
                       sigma_s_d,
                       sigma_f_d);
    
    
}

shared_ptr<Boundary_Source> VERA_Solid_Geometry::
boundary_source(vector<double> const &position) const
{
    Assert(false);
}
