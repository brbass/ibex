#include "Weight_Function.hh"

#include <cmath>
#include <limits>

#include "Angular_Discretization.hh"
#include "Basis_Function.hh"
#include "Cartesian_Plane.hh"
#include "Cross_Section.hh"
#include "Energy_Discretization.hh"
#include "Material.hh"
#include "Meshless_Function.hh"
#include "Quadrature_Rule.hh"
#include "Solid_Geometry.hh"
#include "XML_Node.hh"

using namespace std;
// using std::function;
// using std::max;
// using std::min;
// using std::numeric_limits;
// using std::shared_ptr;
// using std::vector;
namespace qr = Quadrature_Rule;

Weight_Function::
Weight_Function(int index,
                int dimension,
                int integration_ordinates,
                Material_Options material_options,
                shared_ptr<Meshless_Function> meshless_function,
                vector<shared_ptr<Basis_Function> > basis_functions,
                shared_ptr<Solid_Geometry> solid_geometry,
                vector<shared_ptr<Cartesian_Plane> > boundary_surfaces):
    index_(index),
    dimension_(dimension),
    position_(meshless_function->position()),
    integration_ordinates_(integration_ordinates),
    number_of_basis_functions_(basis_functions.size()),
    number_of_boundary_surfaces_(boundary_surfaces.size()),
    radius_(meshless_function->radius()),
    material_options_(material_options),
    meshless_function_(meshless_function),
    basis_functions_(basis_functions),
    solid_geometry_(solid_geometry),
    boundary_surfaces_(boundary_surfaces)
{
    point_type_ = (number_of_boundary_surfaces_ > 0 ?
                   Weight_Function::Point_Type::BOUNDARY :
                   Weight_Function::Point_Type::INTERNAL);

    // Calculate boundary limits
    double lim = numeric_limits<double>::max();
    min_boundary_limits_.assign(dimension_, -lim);
    max_boundary_limits_.assign(dimension_, lim);
    for (int i = 0; i < number_of_boundary_surfaces_; ++i)
    {
        shared_ptr<Cartesian_Plane> surface = boundary_surface(i);
        int dim_sur = surface->surface_dimension();
        double pos_sur = surface->position();
        double n_sur = surface->normal();
        
        if (n_sur < 0)
        {
            min_boundary_limits_[dim_sur] = max(min_boundary_limits_[dim_sur], pos_sur);
        }
        else
        {
            max_boundary_limits_[dim_sur] = min(max_boundary_limits_[dim_sur], pos_sur);
        }
    }
    
    calculate_material();
    calculate_integrals();
    check_class_invariants();
}

bool Weight_Function::
get_full_quadrature(vector<vector<double> > &ordinates,
                    vector<double> &weights) const
{
    switch (dimension_)
    {
    case 1:
        return get_full_quadrature_1d(ordinates,
                                      weights);
    case 2:
        return get_full_quadrature_2d(ordinates,
                                      weights);
    default:
        return false;
    }
}

bool Weight_Function::
get_full_quadrature_1d(vector<vector<double> > &ordinates,
                       vector<double> &weights) const
{
    // Get min and max of weight
    double position = position_[0];
    double x1 = position - radius_;
    double x2 = position + radius_;

    // Compare to boundaries
    x1 = max(x1, min_boundary_limits_[0]);
    x2 = min(x2, max_boundary_limits_[0]);
    
    // Get quadrature
    vector<double> ordinates_x;
    bool success =  qr::cartesian_1d(qr::Quadrature_Type::GAUSS_LEGENDRE,
                                     integration_ordinates_,
                                     x1,
                                     x2,
                                     ordinates_x,
                                     weights);
    qr::convert_to_position_1d(ordinates_x,
                               ordinates);
    return success;
}
    
bool Weight_Function::
get_full_quadrature_2d(vector<vector<double> > &ordinates,
                       vector<double> &weights) const
{
    // Return standard cylindrical quadrature if no boundary surfaces
    bool success = false;
    vector<double> ordinates_x;
    vector<double> ordinates_y;
    if (number_of_boundary_surfaces_ == 0)
    {
        success = qr::cylindrical_2d(qr::Quadrature_Type::GAUSS_LEGENDRE,
                                     qr::Quadrature_Type::GAUSS_LEGENDRE,
                                     integration_ordinates_,
                                     integration_ordinates_,
                                     position_[0],
                                     position_[1],
                                     0.,
                                     radius_,
                                     0.,
                                     2. * M_PI,
                                     ordinates_x,
                                     ordinates_y,
                                     weights);
    }
    else
    {
        success = qr::cartesian_bounded_cylindrical_2d(qr::Quadrature_Type::GAUSS_LEGENDRE,
                                                       qr::Quadrature_Type::GAUSS_LEGENDRE,
                                                       integration_ordinates_,
                                                       integration_ordinates_,
                                                       position_[0],
                                                       position_[1],
                                                       radius_,
                                                       min_boundary_limits_[0],
                                                       max_boundary_limits_[0],
                                                       min_boundary_limits_[1],
                                                       max_boundary_limits_[1],
                                                       ordinates_x,
                                                       ordinates_y,
                                                       weights);
    }

    qr::convert_to_position_2d(ordinates_x,
                               ordinates_y,
                               ordinates);
    
    return success;
}

bool Weight_Function::
get_basis_quadrature(int i,
                     vector<vector<double> > &ordinates,
                     vector<double> &weights) const
{
    switch (dimension_)
    {
    case 1:
        return get_basis_quadrature_1d(i,
                                       ordinates,
                                       weights);
    case 2:
        return get_basis_quadrature_2d(i,
                                       ordinates,
                                       weights);
    default:
        return false;
    }
}

bool Weight_Function::
get_basis_quadrature_1d(int i,
                        vector<vector<double> > &ordinates,
                        vector<double> &weights) const
{
    // Get basis and weight positional information
    shared_ptr<Basis_Function> basis = basis_function(i);
    double basis_position = basis->position()[0];
    double weight_position = position_[0];
    double basis_radius = basis->radius();

    // Get starting/ending positions for basis and weight
    double x1 = max(weight_position - radius_, basis_position - basis_radius);
    double x2 = min(weight_position + radius_, basis_position + basis_radius);

    // Apply boundary surfaces
    if (min_boundary_limits_[0] > x1)
    {
        x1 = min_boundary_limits_[0];
    }
    if (max_boundary_limits_[0] < x2)
    {
        x2 = max_boundary_limits_[0];
    }
    
    // Get quadrature
    vector<double> ordinates_x;
    bool success = qr::cartesian_1d(qr::Quadrature_Type::GAUSS_LEGENDRE,
                                    integration_ordinates_,
                                    x1,
                                    x2,
                                    ordinates_x,
                                    weights);
    qr::convert_to_position_1d(ordinates_x,
                               ordinates);
    return success;
}
    
bool Weight_Function::
get_basis_quadrature_2d(int i,
                        vector<vector<double> > &ordinates,
                        vector<double> &weights) const
{
    // Get basis information
    shared_ptr<Basis_Function> basis = basis_function(i);
    vector<double> basis_position = basis->position();

    // If either weight or basis does not intersect surface,
    // return standard double_cylindrical quadrature
    bool success;
    vector<double> ordinates_x;
    vector<double> ordinates_y;
    if (number_of_boundary_surfaces_ == 0
        || basis->number_of_boundary_surfaces() == 0)
    {
        success = qr::double_cylindrical_2d(qr::Quadrature_Type::GAUSS_LEGENDRE,
                                            qr::Quadrature_Type::GAUSS_LEGENDRE,
                                            integration_ordinates_,
                                            integration_ordinates_,
                                            position_[0],
                                            position_[1],
                                            radius_,
                                            basis_position[0],
                                            basis_position[1],
                                            basis->radius(),
                                            ordinates_x,
                                            ordinates_y,
                                            weights);
    }
    else
    {
        success = qr::cartesian_bounded_double_cylindrical_2d(qr::Quadrature_Type::GAUSS_LEGENDRE,
                                                              qr::Quadrature_Type::GAUSS_LEGENDRE,
                                                              integration_ordinates_,
                                                              integration_ordinates_,
                                                              position_[0],
                                                              position_[1],
                                                              radius_,
                                                              basis_position[0],
                                                              basis_position[1],
                                                              basis->radius(),
                                                              min_boundary_limits_[0],
                                                              max_boundary_limits_[0],
                                                              min_boundary_limits_[1],
                                                              max_boundary_limits_[1],
                                                              ordinates_x,
                                                              ordinates_y,
                                                              weights);
    }

    qr::convert_to_position_2d(ordinates_x,
                               ordinates_y,
                               ordinates);
    
    return success;
}

bool Weight_Function::
get_full_surface_quadrature(int s,
                            vector<vector<double> > &ordinates,
                            vector<double> &weights) const
{
    switch (dimension_)
    {
    case 1:
        ordinates.assign(1, vector<double>({boundary_surfaces_[s]->position()}));
        weights.assign(1, 0.);
        return true;
    case 2:
        return get_full_surface_quadrature_2d(s,
                                              ordinates,
                                              weights);
    default:
        return false;
    }
}

bool Weight_Function::
get_full_surface_quadrature_2d(int s,
                               vector<vector<double> > &ordinates,
                               vector<double> &weights) const
{
    // Get surface information
    shared_ptr<Cartesian_Plane> surface = boundary_surfaces_[s];
    int dim_sur = surface->surface_dimension();
    int dim_other = (dim_sur == 0 ? 1 : 0);
    double pos_sur = surface->position();
    
    // Calculate bounds on surface integral
    double dist = abs(pos_sur - position_[dim_sur]);
    double l = sqrt(radius_ * radius_ - dist * dist);
    double smin = position_[dim_other] - l;
    double smax = position_[dim_other] + l;
    
    // Check for boundaries
    smin = max(smin, min_boundary_limits_[dim_other]);
    smax = min(smax, max_boundary_limits_[dim_other]);

    vector<double> ordinates_main;
    bool success = qr::cartesian_1d(qr::Quadrature_Type::GAUSS_LEGENDRE,
                                    integration_ordinates_,
                                    smin,
                                    smax,
                                    ordinates_main,
                                    weights);
    
    vector<double> ordinates_other(ordinates.size(), pos_sur);
    
    switch (dim_sur)
    {
    case 0:
        qr::convert_to_position_2d(ordinates_other,
                                   ordinates_main,
                                   ordinates);
        return success;
    case 1:
        qr::convert_to_position_2d(ordinates_main,
                                   ordinates_other,
                                   ordinates);
        return success;
    default:
        return false;
    }
}

bool Weight_Function::
get_basis_surface_quadrature(int i,
                             int s,
                             vector<vector<double> > &ordinates,
                             vector<double> &weights) const
{
    switch (dimension_)
    {
    case 1:
        ordinates.assign(1, vector<double>({boundary_surfaces_[s]->position()}));
        weights.assign(1, 0.);
        return true;
    case 2:
        return get_basis_surface_quadrature_2d(i,
                                               s,
                                               ordinates,
                                               weights);
    default:
        return false;
    }
}

bool Weight_Function::
get_basis_surface_quadrature_2d(int i,
                                int s,
                                vector<vector<double> > &ordinates,
                                vector<double> &weights) const
{
    // Get surface information
    shared_ptr<Cartesian_Plane> surface = boundary_surfaces_[s];
    int dim_sur = surface->surface_dimension();
    int dim_other = (dim_sur == 0 ? 1 : 0);
    double pos_sur = surface->position();
    
    // Check basis function boundaries
    vector<double> basis_position = basis_functions_[i]->position();
    double basis_radius = basis_functions_[i]->radius();
    int number_of_basis_boundary_surfaces = basis_functions_[i]->number_of_boundary_surfaces();
    double distb = abs(pos_sur - basis_position[dim_sur]);
    
    // If basis function does not intersect, return empty quadrature
    if (number_of_basis_boundary_surfaces == 0 || distb > basis_radius)
    {
        ordinates.resize(0);
        weights.resize(0);
        
        return true;
    }

    double lb = sqrt(basis_radius * basis_radius - distb * distb);
    double sbmin = basis_position[dim_other] - lb;
    double sbmax = basis_position[dim_other] + lb;
    
    // Calculate bounds on surface integral
    double dist = abs(pos_sur - position_[dim_sur]);
    double l = sqrt(radius_ * radius_ - dist * dist);
    double smin = max(position_[dim_other] - l, sbmin);
    double smax = min(position_[dim_other] + l, sbmax);
    
    // Check for boundaries
    smin = max(smin, min_boundary_limits_[dim_other]);
    smax = min(smax, max_boundary_limits_[dim_other]);

    if (smin > smax)
    {
        ordinates.resize(0);
        weights.resize(0);

        return true;
    }
    
    vector<double> ordinates_main;
    bool success = qr::cartesian_1d(qr::Quadrature_Type::GAUSS_LEGENDRE,
                                    integration_ordinates_,
                                    smin,
                                    smax,
                                    ordinates_main,
                                    weights);
    
    vector<double> ordinates_other(ordinates.size(), pos_sur);
    
    switch (dim_sur)
    {
    case 0:
        qr::convert_to_position_2d(ordinates_other,
                                   ordinates_main,
                                   ordinates);
        return success;
    case 1:
        qr::convert_to_position_2d(ordinates_main,
                                   ordinates_other,
                                   ordinates);
        return success;
    default:
        return false;
    }
}

void Weight_Function::
calculate_integrals()
{
    // Perform basis/weight integrals
    is_b_w_.assign(number_of_boundary_surfaces_ * number_of_basis_functions_, 0);
    iv_b_w_.assign(number_of_basis_functions_, 0);
    iv_b_dw_.assign(number_of_basis_functions_ * dimension_, 0);
    iv_db_w_.assign(number_of_basis_functions_ * dimension_, 0);
    iv_db_dw_.assign(number_of_basis_functions_ * dimension_ * dimension_, 0);

    shared_ptr<Meshless_Function> weight = meshless_function_;
    
    for (int i = 0; i < number_of_basis_functions_; ++i)
    {
        shared_ptr<Meshless_Function> basis = basis_functions_[i]->function();

        // Surface integrals
        for (int s = 0; s < number_of_boundary_surfaces_; ++s)
        {
            // Get quadrature
            vector<vector<double> > ordinates;
            vector<double> weights;
            Assert(get_basis_surface_quadrature(i,
                                                s,
                                                ordinates,
                                                weights));
            int number_of_ordinates = ordinates.size();
            
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                vector<double> position = ordinates[o];
                double b = basis->value(position);
                double w = weight->value(position);

                is_b_w_[s + number_of_boundary_surfaces_ * i] += weights[o] * b * w;
            }
        }

        // Volume integrals
        vector<vector<double> > ordinates;
        vector<double> weights;
        
        Assert(get_basis_quadrature(i,
                                    ordinates,
                                    weights));
        int number_of_ordinates = ordinates.size();

        for (int o = 0; o < number_of_ordinates; ++o)
        {
            vector<double> position = ordinates[o];
            double b = basis->value(position);
            double w = weight->value(position);
            vector<double> db = basis->gradient_value(position);
            vector<double> dw = weight->gradient_value(position);
            
            iv_b_w_[i] += weights[o] * b * w;

            for (int d1 = 0; d1 < dimension_; ++d1)
            {
                int k1 = d1 + dimension_ * i;
                iv_b_dw_[k1] += weights[o] * b * dw[d1];
                iv_db_w_[k1] += weights[o] * b * dw[d1];
                
                for (int d2 = 0; d2 < dimension_; ++d2)
                {
                    int k2 = d1 + dimension_ * (d2 + dimension_ * i);
                    iv_db_dw_[k2] += weights[o] * db[d1] * dw[d2];
                }
            }
        }
    }
}

void Weight_Function::
calculate_material()
{
    // Make sure options for material weighting are cohesive
    switch (material_options_.weighting)
    {
    case Material_Options::Weighting::POINT:
        material_options_.total = Material_Options::Total::ISOTROPIC;
        break; // POINT
    case Material_Options::Weighting::WEIGHT:
        material_options_.total = Material_Options::Total::ISOTROPIC;
        break; // WEIGHT
    case Material_Options::Weighting::FLUX:
        AssertMsg(false, "Material_Options.flux must be set in order to weight with flux");
        break; // FLUX
    }
    
    switch (material_options_.output)
    {
    case Material_Options::Output::STANDARD:
    {
        switch(material_options_.weighting)
        {
        case Material_Options::Weighting::POINT:
        {
            return calculate_standard_point_material();
        } // POINT
        case Material_Options::Weighting::WEIGHT:
        {
            return calculate_standard_weight_material();
        } // WEIGHT
        case Material_Options::Weighting::FLUX:
        {
            AssertMsg(false, "Weight_Function flux weighting not yet implemented");
            return;
        } // FLUX
        } // Weighting
        break;
    } // STANDARD
    case Material_Options::Output::SUPG:
    {
        switch(material_options_.weighting)
        {
        case Material_Options::Weighting::POINT:
        {
            return calculate_supg_point_material();
        } // POINT
        case Material_Options::Weighting::WEIGHT:
        {
            return calculate_supg_weight_material();
        } // WEIGHT
        case Material_Options::Weighting::FLUX:
        {
            AssertMsg(false, "Weight_Function flux weighting not yet implemented");
            return;
        } // FLUX
        } // Weighting
        break;
    } // SUPG
    } // Output
}

void Weight_Function::
calculate_standard_point_material()
{
    // Get ordinates
    vector<vector<double> > ordinates;
    vector<double> weights;
    Assert(get_full_quadrature(ordinates,
                               weights));
    int number_of_ordinates = ordinates.size();
    
    // Get material for angular and energy discretization
    shared_ptr<Material> test_material = solid_geometry_->material(position_);
    shared_ptr<Angular_Discretization> angular_discretization = test_material->angular_discretization();
    shared_ptr<Energy_Discretization> energy_discretization = test_material->energy_discretization();
    
    // Get material data
    int number_of_groups = energy_discretization->number_of_groups();
    int number_of_scattering_moments = angular_discretization->number_of_scattering_moments();
    int number_of_moments = angular_discretization->number_of_moments();
    
    // Calculate weight function integral
    double weight_integral = 0.;
    for (int i = 0; i < number_of_ordinates; ++i)
    {
        vector<double> position = ordinates[i];
        double w = meshless_function_->value(position);
        weight_integral += w * weights[i];
    }

    // Get the weighted internal source
    vector<double> const internal_source_temp = test_material->internal_source()->data();
    vector<double> internal_source_v(number_of_groups, 0.);
    for (int g = 0; g < number_of_groups; ++g)
    {
        internal_source_v[g] = weight_integral * internal_source_temp[g];
    }
    
    shared_ptr<Cross_Section> internal_source
        = make_shared<Cross_Section>(test_material->internal_source()->dependencies(),
                                     angular_discretization,
                                     energy_discretization,
                                     internal_source_v);
    
    material_
        = make_shared<Material>(index_,
                                angular_discretization,
                                energy_discretization,
                                test_material->sigma_t(),
                                test_material->sigma_s(),
                                test_material->nu(),
                                test_material->sigma_f(),
                                test_material->chi(),
                                internal_source);
}

void Weight_Function::
calculate_standard_weight_material()
{
    // Get ordinates
    vector<vector<double> > ordinates;
    vector<double> weights;
    Assert(get_full_quadrature(ordinates,
                               weights));
    int number_of_ordinates = ordinates.size();
    
    // Get material for angular and energy discretization
    shared_ptr<Material> test_material = solid_geometry_->material(position_);
    shared_ptr<Angular_Discretization> angular_discretization = test_material->angular_discretization();
    shared_ptr<Energy_Discretization> energy_discretization = test_material->energy_discretization();
    
    // Get material data
    int number_of_groups = energy_discretization->number_of_groups();
    int number_of_scattering_moments = angular_discretization->number_of_scattering_moments();
    int number_of_moments = angular_discretization->number_of_moments();
    
    vector<double> sigma_t_v(number_of_groups, 0.);
    vector<double> sigma_s_v(number_of_groups * number_of_groups * number_of_scattering_moments, 0.);
    vector<double> nu_v(number_of_groups, 0.);
    vector<double> sigma_f_v(number_of_groups, 0.);
    vector<double> chi_v(number_of_groups, 0.);
    vector<double> internal_source_v(number_of_groups, 0.);
    double weight_integral = 0.;

    // Calculate numerator and denominator separately
    for (int i = 0; i < number_of_ordinates; ++i)
    {
        vector<double> position = ordinates[i];
        shared_ptr<Material> material = solid_geometry_->material(position);
        vector<double> const sigma_t_temp = material->sigma_t()->data();
        vector<double> const sigma_s_temp = material->sigma_t()->data();
        vector<double> const nu_temp = material->sigma_t()->data();
        vector<double> const sigma_f_temp = material->sigma_t()->data();
        vector<double> const chi_temp = material->sigma_t()->data();
        vector<double> const internal_source_temp = material->internal_source()->data();
        double w = meshless_function_->value(position);
        weight_integral += w * weights[i];
                
        for (int g = 0; g < number_of_groups; ++g)
        {
            sigma_t_v[g] += w * sigma_t_temp[g] * weights[i];
            sigma_s_v[g] += w * sigma_s_temp[g] * weights[i];
            nu_v[g] += w * nu_temp[g] * weights[i];
            sigma_f_v[g] += w * sigma_f_temp[g] * weights[i];
            chi_v[g] += w * chi_temp[g] * weights[i];
            internal_source_v[g] += w * internal_source_temp[g] * weights[i];
                    
            for (int g2 = 0; g2 < number_of_groups; ++g2)
            {
                for (int m = 0; m < number_of_scattering_moments; ++m)
                {
                    int k = g + number_of_groups * (g2 + number_of_groups * m);
                    sigma_s_v[k] += w * sigma_s_temp[k] * weights[i];
                }
            }
        }
    }
            
    // Divide numerator by denominator for cross sections
    for (int g = 0; g < number_of_groups; ++g)
    {
        sigma_t_v[g] /= weight_integral;
        sigma_s_v[g] /= weight_integral;
        nu_v[g] /= weight_integral;
        sigma_f_v[g] /= weight_integral;
        chi_v[g] /= weight_integral;
                
        for (int g2 = 0; g2 < number_of_groups; ++g2)
        {
            for (int m = 0; m < number_of_scattering_moments; ++m)
            {
                int k = g + number_of_groups * (g2 + number_of_groups * m);
                sigma_s_v[k] /= weight_integral;
            }
        }
    }

    // Create cross sections
    Cross_Section::Dependencies none_group;
    none_group.energy = Cross_Section::Dependencies::Energy::GROUP;
    Cross_Section::Dependencies scattering_group2;
    scattering_group2.angular = Cross_Section::Dependencies::Angular::SCATTERING_MOMENTS;
    scattering_group2.energy = Cross_Section::Dependencies::Energy::GROUP_TO_GROUP;
    
    shared_ptr<Cross_Section> sigma_t
        = make_shared<Cross_Section>(none_group,
                                     angular_discretization,
                                     energy_discretization,
                                     sigma_t_v);
    shared_ptr<Cross_Section> sigma_s
        = make_shared<Cross_Section>(scattering_group2,
                                     angular_discretization,
                                     energy_discretization,
                                     sigma_s_v);
    shared_ptr<Cross_Section> nu
        = make_shared<Cross_Section>(none_group,
                                     angular_discretization,
                                     energy_discretization,
                                     nu_v);
    shared_ptr<Cross_Section> sigma_f
        = make_shared<Cross_Section>(none_group,
                                     angular_discretization,
                                     energy_discretization,
                                     sigma_f_v);
    shared_ptr<Cross_Section> chi
        = make_shared<Cross_Section>(none_group,
                                     angular_discretization,
                                     energy_discretization,
                                     chi_v);
    shared_ptr<Cross_Section> internal_source
        = make_shared<Cross_Section>(none_group,
                                     angular_discretization,
                                     energy_discretization,
                                     internal_source_v);
    
    // Create material
    material_ = make_shared<Material>(index_,
                                      angular_discretization,
                                      energy_discretization,
                                      sigma_t,
                                      sigma_s,
                                      nu,
                                      sigma_f,
                                      chi,
                                      internal_source);
    
}

void Weight_Function::
calculate_supg_point_material()
{
    // Get ordinates
    vector<vector<double> > ordinates;
    vector<double> weights;
    Assert(get_full_quadrature(ordinates,
                               weights));
    int number_of_ordinates = ordinates.size();
    
    // Get material for angular and energy discretization
    shared_ptr<Material> test_material = solid_geometry_->material(position_);
    shared_ptr<Angular_Discretization> angular_discretization = test_material->angular_discretization();
    shared_ptr<Energy_Discretization> energy_discretization = test_material->energy_discretization();
    
    // Get material data
    int number_of_groups = energy_discretization->number_of_groups();
    int number_of_scattering_moments = angular_discretization->number_of_scattering_moments();
    int number_of_moments = angular_discretization->number_of_moments();
    int dimensionp1 = dimension_ + 1;
    
    // Calculate weight function integral
    vector<double> weight_integral(dimension_ + 1, 0.);
    for (int i = 0; i < number_of_ordinates; ++i)
    {
        vector<double> position = ordinates[i];
        
        // Get function weight values
        double w1 = meshless_function_->value(position);
        vector<double> w2 = meshless_function_->gradient_value(position);
        vector<double> w(dimension_ + 1, 0);
        w[0] = w1;
        for (int j = 0; j < dimension_; ++j)
        {
            w[j+1] = w2[j];
        }
                
        for (int j = 0; j < dimension_ + 1; ++j)
        {
            weight_integral[j] += w[j] * weights[i];
        }
    }
            
    // Get the weighted internal source
    vector<double> const internal_source_temp = test_material->internal_source()->data();
    vector<double> internal_source_v(dimensionp1 * number_of_groups, 0.);
    for (int j = 0; j < dimension_ + 1; ++j)
    {    
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k = j + dimensionp1 * g;
            internal_source_v[k] = internal_source_temp[g] * weight_integral[j];
        }
    }

    Cross_Section::Dependencies dep = test_material->internal_source()->dependencies();
    dep.dimensional = Cross_Section::Dependencies::Dimensional::SUPG;
    shared_ptr<Cross_Section> internal_source
        = make_shared<Cross_Section>(dep,
                                     angular_discretization,
                                     energy_discretization,
                                     internal_source_v);
    
    // Create material
    material_ = make_shared<Material>(index_,
                                      angular_discretization,
                                      energy_discretization,
                                      test_material->sigma_t(),
                                      test_material->sigma_s(),
                                      test_material->nu(),
                                      test_material->sigma_f(),
                                      test_material->chi(),
                                      internal_source);
}

void Weight_Function::
calculate_supg_weight_material()
{
    // Get ordinates
    vector<vector<double> > ordinates;
    vector<double> weights;
    Assert(get_full_quadrature(ordinates,
                               weights));
    int number_of_ordinates = ordinates.size();
    
    // Get material for angular and energy discretization
    shared_ptr<Material> test_material = solid_geometry_->material(position_);
    shared_ptr<Angular_Discretization> angular_discretization = test_material->angular_discretization();
    shared_ptr<Energy_Discretization> energy_discretization = test_material->energy_discretization();
    
    // Get material data
    int number_of_groups = energy_discretization->number_of_groups();
    int number_of_scattering_moments = angular_discretization->number_of_scattering_moments();
    int number_of_moments = angular_discretization->number_of_moments();
    int dimensionp1 = dimension_ + 1;
    
    vector<double> sigma_t_v(dimensionp1 * number_of_groups, 0);
    vector<double> sigma_s_v(dimensionp1 * number_of_groups * number_of_groups * number_of_scattering_moments, 0);
    vector<double> nu_v(dimensionp1 * number_of_groups, 0);
    vector<double> sigma_f_v(dimensionp1 * number_of_groups, 0);
    vector<double> chi_v(dimensionp1 * number_of_groups, 0);
    vector<double> internal_source_v(dimensionp1 * number_of_groups, 0);
    vector<double> weight_integral(dimensionp1, 0.);
            
    // Calculate numerator and denominator separately
    for (int i = 0; i < number_of_ordinates; ++i)
    {
        vector<double> position = ordinates[i];
        shared_ptr<Material> material = solid_geometry_->material(position);
        vector<double> const sigma_t_temp = material->sigma_t()->data();
        vector<double> const sigma_s_temp = material->sigma_t()->data();
        vector<double> const nu_temp = material->sigma_t()->data();
        vector<double> const sigma_f_temp = material->sigma_t()->data();
        vector<double> const chi_temp = material->sigma_t()->data();
        vector<double> const internal_source_temp = material->internal_source()->data();

        // Get function weight values
        double w1 = meshless_function_->value(position);
        vector<double> w2 = meshless_function_->gradient_value(position);
        vector<double> w(dimension_ + 1, 0);
        w[0] = w1;
        for (int j = 0; j < dimension_; ++j)
        {
            w[j+1] = w2[j];
        }
                
        for (int j = 0; j < dimension_ + 1; ++j)
        {
            weight_integral[j] += w[j] * weights[i];
            
            for (int g = 0; g < number_of_groups; ++g)
            {
                int k1 = j + dimensionp1 * g;
                sigma_t_v[k1] += w[j] * sigma_t_temp[g] * weights[i];
                sigma_s_v[k1] += w[j] * sigma_s_temp[g] * weights[i];
                nu_v[k1] += w[j] * nu_temp[g] * weights[i];
                sigma_f_v[k1] += w[j] * sigma_f_temp[g] * weights[i];
                chi_v[k1] += w[j] * chi_temp[g] * weights[i];
                internal_source_v[k1] += w[j] * internal_source_temp[g] * weights[i];
                        
                for (int g2 = 0; g2 < number_of_groups; ++g2)
                {
                    for (int m = 0; m < number_of_scattering_moments; ++m)
                    {
                        int k2 = j + dimensionp1 * (g + number_of_groups * (g2 + number_of_groups * m));
                        int k3 = g + number_of_groups * (g2 + number_of_groups * m);
                        sigma_s_v[k2] += w[j] * sigma_s_temp[k3] * weights[i];
                    }
                }
            }
        }
    }
            
    // Divide numerator by denominator for cross sections
    for (int j = 0; j < dimension_ + 1; ++j)
    {    
        for (int g = 0; g < number_of_groups; ++g)
        {
            int k1 = j + dimensionp1 * g;
            sigma_t_v[k1] /= weight_integral[j];
            sigma_s_v[k1] /= weight_integral[j];
            nu_v[k1] /= weight_integral[j];
            sigma_f_v[k1] /= weight_integral[j];
            chi_v[k1] /= weight_integral[j];
            
            for (int g2 = 0; g2 < number_of_groups; ++g2)
            {
                for (int m = 0; m < number_of_scattering_moments; ++m)
                {
                    int k2 = j + dimensionp1 * (g + number_of_groups * (g2 + number_of_groups * m));
                    sigma_s_v[k2] /= weight_integral[j];
                }
            }
        }
    }
            
    // Create cross sections
    Cross_Section::Dependencies none_group;
    none_group.energy = Cross_Section::Dependencies::Energy::GROUP;
    none_group.dimensional = Cross_Section::Dependencies::Dimensional::SUPG;
    Cross_Section::Dependencies scattering_group2;
    scattering_group2.angular = Cross_Section::Dependencies::Angular::SCATTERING_MOMENTS;
    scattering_group2.energy = Cross_Section::Dependencies::Energy::GROUP_TO_GROUP;
    scattering_group2.dimensional = Cross_Section::Dependencies::Dimensional::SUPG;

    shared_ptr<Cross_Section> sigma_t
        = make_shared<Cross_Section>(none_group,
                                     angular_discretization,
                                     energy_discretization,
                                     sigma_t_v);
    shared_ptr<Cross_Section> sigma_s
        = make_shared<Cross_Section>(scattering_group2,
                                     angular_discretization,
                                     energy_discretization,
                                     sigma_s_v);
    shared_ptr<Cross_Section> nu
        = make_shared<Cross_Section>(none_group,
                                     angular_discretization,
                                     energy_discretization,
                                     nu_v);
    shared_ptr<Cross_Section> sigma_f
        = make_shared<Cross_Section>(none_group,
                                     angular_discretization,
                                     energy_discretization,
                                     sigma_f_v);
    shared_ptr<Cross_Section> chi
        = make_shared<Cross_Section>(none_group,
                                     angular_discretization,
                                     energy_discretization,
                                     chi_v);
    shared_ptr<Cross_Section> internal_source
        = make_shared<Cross_Section>(none_group,
                                     angular_discretization,
                                     energy_discretization,
                                     internal_source_v);
    
    // Create material
    material_ = make_shared<Material>(index_,
                                      angular_discretization,
                                      energy_discretization,
                                      sigma_t,
                                      sigma_s,
                                      nu,
                                      sigma_f,
                                      chi,
                                      internal_source);
}

void Weight_Function::
calculate_boundary_source()
{
    switch(material_options_.weighting)
    {
    case Material_Options::Weighting::POINT:
        return calculate_point_boundary_source();
    case Material_Options::Weighting::WEIGHT:
        return calculate_weight_boundary_source();
    case Material_Options::Weighting::FLUX:
        AssertMsg(false, "Weight_Function flux weighting not yet implemented");
        return;
    }
}

void Weight_Function::
calculate_point_boundary_source()
{
    AssertMsg(false, "Not yet implemented");
}

void Weight_Function::
calculate_weight_boundary_source()
{
    AssertMsg(false, "Not yet implemented");
}

void Weight_Function::
check_class_invariants() const
{
    // Check that all boundary surfaces actually intersect
    
}

void Weight_Function::
output(XML_Node output_node) const
{
    
}
