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
get_full_quadrature_1d(vector<double> &ordinates,
                       vector<double> &weights) const
{
    // Get min and max of weight
    double position = position_[0];
    double x1 = position - radius_;
    double x2 = position + radius_;

    // Compare to boundaries
    if (x1 < boundary_limits[0]
    for (int i = 0; i < number_of_boundary_surfaces_; ++i)
    {
        shared_ptr<Cartesian_Plane> surface = boundary_surface(i);
        double x_sur = surface->position();
        double n_sur = surface->normal();
        
        if (n_sur < 0)
        {
            if (x_sur > x1)
            {
                x1 = x_sur;
            }
        }
        else
        {
            if (x_sur < x2)
            {
                x2 = x_sur;
            }
        }
    }

    // Get quadrature
    return qr::cartesian_1d(qr::Quadrature_Type::GAUSS_LEGENDRE,
                            integration_ordinates_,
                            x1,
                            x2,
                            ordinates,
                            weights);
}

bool Weight_Function::
get_full_quadrature_2d(vector<double> &ordinates_x,
                       vector<double> &ordinates_y,
                       vector<double> &weights) const
{
    // Return standard cylindrical quadrature if no boundary surfaces
    if (number_of_boundary_surfaces_ == 0)
    {
        return qr::cylindrical_2d(qr::Quadrature_Type::GAUSS_LEGENDRE,
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

    // Get quadrature
    return qr::cartesian_bounded_cylindrical_2d(qr::Quadrature_Type::GAUSS_LEGENDRE,
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

bool Weight_Function::
get_basis_quadrature_1d(int i,
                        vector<double> &ordinates,
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
    return qr::cartesian_1d(qr::Quadrature_Type::GAUSS_LEGENDRE,
                            integration_ordinates_,
                            x1,
                            x2,
                            ordinates,
                            weights);
}
    
bool Weight_Function::
get_basis_quadrature_2d(int i,
                        vector<double> &ordinates_x,
                        vector<double> &ordinates_y,
                        vector<double> &weights) const
{
    // Get basis information
    shared_ptr<Basis_Function> basis = basis_function(i);
    vector<double> basis_position = basis->position();

    // If either weight or basis does not intersect surface,
    // return standard double_cylindrical quadrature
    if (number_of_boundary_surfaces_ == 0
        || basis->number_of_boundary_surfaces() == 0)
    {
        return qr::double_cylindrical_2d(qr::Quadrature_Type::GAUSS_LEGENDRE,
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

    return qr::cartesian_bounded_double_cylindrical_2d(qr::Quadrature_Type::GAUSS_LEGENDRE,
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

void Weight_Function::
calculate_integrals_1d()
{
    // Initialize integrals
    is_b_w_.assign(number_of_boundary_surfaces_ * number_of_basis_functions_, 0);
    iv_b_w_.assign(number_of_basis_functions_, 0);
    iv_b_dw_.assign(number_of_basis_functions_ * dimension_, 0);
    iv_db_dw_.assign(number_of_basis_functions_ * dimension_ * dimension_, 0);
    for (int i = 0; i < number_of_basis_functions_; ++i)
    {
        shared_ptr<Meshless_Function> basis = basis_functions_[i]->function();
        
        // Surface integrals
        for (int s = 0; s < number_of_boundary_surfaces_; ++s)
        {
            vector<double> position = {boundary_surfaces_[i]->position()};
            double b = basis->basis(position);
            double w = meshless_function_->basis(position);

            is_b_w_[s + number_of_boundary_surfaces_ * i] = b * w;
        }
        
        // Volume integrals
        vector<double> ordinates;
        vector<double> weights;
        
        Assert(get_basis_quadrature_1d(i,
                                       ordinates,
                                       weights));
        int number_of_ordinates = ordinates.size();
        
        for (int o = 0; o < number_of_ordinates; ++o)
        {
            vector<double> position = {ordinates[o]};
            double b = basis->basis(position);
            double w = meshless_function_->basis(position);
            double db = basis->d_basis(0, // dim
                                       position);
            double dw = meshless_function_->d_basis(0, // dim
                                                    position);
            
            iv_b_w_[i] += weights[o] * b * w;
            iv_b_dw[i] += weights[o] * b * dw;
            iv_db_dw[i] += weights[o] * db * dw;
        }
    }
}

void Weight_Function::
calculate_integrals_2d()
{
    // Initialize integrals
    is_b_w_.assign(number_of_boundary_surfaces_ * number_of_basis_functions_, 0);
    iv_b_w_.assign(number_of_basis_functions_, 0);
    iv_b_dw_.assign(number_of_basis_functions_ * dimension_, 0);
    iv_db_dw_.assign(number_of_basis_functions_ * dimension_ * dimension_, 0);
    for (int i = 0; i < number_of_basis_functions_; ++i)
    {
        shared_ptr<Meshless_Function> basis = basis_functions_[i]->function();
        vector<double> basis_position = basis_functions_[i]->position();
        double basis_radius = basis_functions_[i]->radius();
        int number_of_basis_boundary_surfaces = basis_functions_[i]->number_of_boundary_surfaces();

        if (number_of_basis_boundary_surfaces != 0)
        {
            // Surface integrals
            for (int s = 0; s < number_of_boundary_surfaces_; ++s)
            {
                shared_ptr<Cartesian_Plane> surface = boundary_surfaces_[s1];
                int dim_sur = surface->surface_dimension();
                double pos_sur = surface->position();
                double n_sur = surface->normal();
                
                // Check if basis function intersects this surface
                if (abs(basis_position[dim_sur] - pos_sur) < basis_radius)
                {
                    // Calculate range of surface integral
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
