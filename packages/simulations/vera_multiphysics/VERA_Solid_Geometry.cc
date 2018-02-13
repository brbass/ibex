#include "VERA_Solid_Geometry.hh"

#include <cmath>

#include "Angular_Discretization.hh"
#include "Energy_Discretization.hh"
#include "Boundary_Source.hh"
#include "Cartesian_Plane.hh"
#include "Cross_Section.hh"
#include "Material.hh"
#include "Material_Factory.hh"

using namespace std;

VERA_Solid_Geometry::
VERA_Solid_Geometry(bool include_ifba,
                    std::shared_ptr<VERA_Temperature> temperature,
                    shared_ptr<Angular_Discretization> angular,
                    shared_ptr<Energy_Discretization> energy,
                    vector<shared_ptr<Material> > materials,
                    shared_ptr<Boundary_Source> boundary_source):
    include_ifba_(include_ifba),
    temperature_(temperature),
    angular_(angular),
    energy_(energy),
    material_factory_(make_shared<Material_Factory>(angular,
                                                    energy)),
    number_of_material_types_(5),
    number_of_problem_types_(2),
    number_of_temperature_types_(2),
    materials_(materials),
    boundary_source_(boundary_source)
{
    // Set up radii
    material_radii2_
        = {0.4096, 0.4106, 0.418, 0.475, 0.891};
    for (double &radius : material_radii2_)
    {
        radius = radius * radius;
    }

    // Get boundary surfaces
    double prob_radius = 0.63;
    boundary_surfaces_.resize(4);
    for (int d = 0; d < 2; ++d)
    {
        for (int n = 0; n < 2; ++n)
        {
            double normal = n == 0 ? -1 : 1;
            double position = n == 0 ? -prob_radius : prob_radius;
            int index = n + 2 * d;
            boundary_surfaces_[index]
                = make_shared<Cartesian_Plane>(index,
                                               2, // dimension
                                               Surface::Surface_Type::BOUNDARY,
                                               d, // surface dimension
                                               position,
                                               normal);
            boundary_surfaces_[index]->set_boundary_source(boundary_source);
        }
    }
    
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
    Assert(angular_);
    Assert(energy_);
    Assert(material_factory_);
    Assert(angular_->dimension() == 2);
    Assert(angular_->number_of_scattering_moments() == 2);
    Assert(energy_->number_of_groups() == 2);
    Assert(materials_.size() == number_of_material_types_ * number_of_problem_types_ * number_of_temperature_types_);
}

shared_ptr<Material> VERA_Solid_Geometry::
get_material_by_index(Material_Type mat_type,
                      Problem_Type prob_type,
                      Temperature_Type temp_type) const
{
    int index = static_cast<int>(mat_type) + number_of_material_types_ * (static_cast<int>(prob_type) + number_of_problem_types_ * static_cast<int>(temp_type));
    
    return materials_[index];
}

double VERA_Solid_Geometry::
radial_distance2(vector<double> const &position) const
{
    return position[0] * position[0] + position[1] * position[1];
}

shared_ptr<Material> VERA_Solid_Geometry::
material(vector<double> const &position) const
{
    // Get material by position
    double const radius2 = radial_distance2(position);
    Material_Type mat_type = NONE;
    for (int i = 0; i < number_of_material_types_; ++i)
    {
        if (radius2 < material_radii2_[i])
        {
            mat_type = static_cast<Material_Type>(i);
            break;
        }
    }
    Assert(mat_type != NONE);

    // Get problem type
    Problem_Type prob_type = include_ifba_ ? V1E : V1B;

    // Get materials at two temperatures
    shared_ptr<Material> mat600 = get_material_by_index(mat_type,
                                                        prob_type,
                                                        K600);
    shared_ptr<Material> mat1200 = get_material_by_index(mat_type,
                                                         prob_type,
                                                         K1200);
    
    // Weight materials according to temperature
    double const temperature = (*temperature_)(position);
    
    return weighted_material(temperature,
                             mat600,
                             mat1200);
}

vector<double> VERA_Solid_Geometry::
weighted_cross_section(double temperature,
                       vector<double> const &xs600,
                       vector<double> const &xs1200) const
{
    double const t600 = 600;
    double const t1200 = 1200;
    double const t600_t = sqrt(t600 / temperature);
    double const t600_1200 = sqrt(t600 / t1200);
    double const mult600 = t600_1200 - t600_t;
    double const mult1200 = t600_t - 1;
    double const multden = 1. / (t600_1200 - 1);
    
    int size = xs600.size();
    vector<double> xsnew(size);
    for (int i = 0; i < size; ++i)
    {
        xsnew[i] = (xs600[i] * mult600 + xs1200[i] * mult1200) * multden;
    }
    
    return xsnew;
}

shared_ptr<Material> VERA_Solid_Geometry::
weighted_material(double temperature,
                  shared_ptr<Material> mat600,
                  shared_ptr<Material> mat1200) const
{
    // Get weighted cross sections
    vector<double> sigma_t
        = weighted_cross_section(temperature,
                                 mat600->sigma_t()->data(),
                                 mat1200->sigma_t()->data());
    vector<double> sigma_s
        = weighted_cross_section(temperature,
                                 mat600->sigma_s()->data(),
                                 mat1200->sigma_s()->data());
    vector<double> sigma_f
        = weighted_cross_section(temperature,
                                 mat600->sigma_f()->data(),
                                 mat1200->sigma_f()->data());
    vector<double> internal_source
        = weighted_cross_section(temperature,
                                 mat600->internal_source()->data(),
                                 mat1200->internal_source()->data());

    // Get material
    return material_factory_->get_full_fission_material(0, // index
                                                        sigma_t,
                                                        sigma_s,
                                                        sigma_f,
                                                        internal_source);
}

shared_ptr<Boundary_Source> VERA_Solid_Geometry::
boundary_source(vector<double> const &position) const
{
    return boundary_source_;
}
