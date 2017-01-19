#include "Weight_Function.hh"

#include <limits>

#include "Basis_Function.hh"
#include "Cartesian_Plane.hh"
#include "Meshless_Function.hh"
#include "Solid_Geometry.hh"
#include "Quadrature_Rule.hh"

using std::numeric_limits;
using std::shared_ptr;
using std::vector;
namspace qr = Quadrature_Rule;

Weight_Function::
Weight_Function(int index,
                int dimension,
                int integration_ordinates,
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
            if (min_boundary_limits[dim_sur] < pos_sur)
            {
                min_boundary_limits[dim_sur] = pos_sur;
            }
        }
        else
        {
            if (max_boundary_limits[dim_sur] > pos_sur)
            {
                max_boundary_limits[dim_sur] = pos_sur;
            }
        }
    }

    check_class_invariants();
    calculate_material();
    calculate_integrals();
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
                                     ordinates_y);
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
        vector<double> position = {boundary_surfaces_[s]->position()};
        ordinates.assign(1, position);
        weights.assign(1, 0.);
        return true;
    case 2:
        return get_full_surface_quadrature_2d(s,
                                              ordinates,
                                              weights);
    default:
        return false;
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
    
    Vector<double> ordinates_other(ordinates.size(), pos_sur);
    
    switch (dim_sur)
    {
    case 0:
        qr::convert_to_position(ordinates_other,
                                ordinates_main,
                                ordinates);
        return success;
    case 1:
        qr::convert_to_position(ordinates_main,
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
        vector<double> position = {boundary_surfaces_[s]->position()};
        ordinates.assign(1, position);
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
    double distb = abs(pos_sur - basis_position_[dim_sur]);
    
    // If basis function does not intersect, return empty quadrature
    if (number_of_basis_boundary_surfaces == 0 || distb > basis_radius)
    {
        ordinates.resize(0);
        weights.resize(0);
        
        return true;
    }

    double lb = sqrt(basis_radius * basis_radius - dist * dist);
    double sbmin = basis_position_[dim_other] - lb;
    double sbmax = basis_position_[dim_other] + lb;
    
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
        qr::convert_to_position(ordinates_other,
                                ordinates_main,
                                ordinates);
        return success;
    case 1:
        qr::convert_to_position(ordinates_main,
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
    // Initialize integrals
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
            vector<double> position = ordiantes[o];
            double b = basis->value(position);
            double w = weight->value(position);
            vector<double> db = basis->gradient_value(position);
            vector<double> dw = weight->gradient_value(position);
            
            iv_b_w_[i] += weights[o] * b * w;

            for (int d1 = 0; d1 < dimension_; ++d1)
            {
                int k1 = d1 + dimension_ * i;
                iv_b_dw[k1] += weights[o] * b * dw[d1];
                iv_db_w[k1] += weights[o] * b * dw[d1];
                
                for (int d2 = 0; d2 < dimension_; ++d2)
                {
                    int k2 = d1 + dimension_ * (d2 + dimension_ * i);
                    iv_db_dw[k2] += weights[o] * db[d1] * dw[d2];
                }
            }
        }
    }
}

void Weight_Function::
check_class_invariants() const
{
    // Check that all boundary surfaces actually intersect
    
}
