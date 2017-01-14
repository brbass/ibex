#include "Weight_Function.hh"

#include "Basis_Function.hh"
#include "Cartesian_Plane.hh"
#include "Meshless_Function.hh"
#include "Solid_Geometry.hh"

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
    integration_ordinates_(integration_ordinates),
    number_of_basis_functions_(basis_functions.size()),
    number_of_boundary_surfaces_(boundary_surfaces.size()),
    
    meshless_function_(meshless_function),
    basis_functions_(basis_functions),
    solid_geometry_(solid_geometry),
    boundary_surfaces_(boundary_surfaces)
{
}

void Weight_Function::
get_full_quadrature_1d(vector<double> &quadrature,
                       vector<double> &weight)
{
    double radius = function()->radius();
    double position = function()->position()[0];
    double x1 = position - radius;
    double x2 = position + radius;
    
    for (int i = 0; i < number_of_boundary_surfaces_; ++i)
    {
        shared_ptr<Cartesian_Plane> surface = boundary_surfaces(i);
        double x_sur = surface->position();
        double n_sur = surface->normal();
        
        if (n_sur > 0)
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
    
    Assert(cartesian_1d(Quadrature_Rule::Quadrature_Type::GAUSS_LEGENDRE,
                        integration_ordinates,
                        x1,
                        x2,
                        ordinates,
                        weights));
}

void Weight_Function::
get_basis_quadrature_1d(int i,
                        vector<double> &quadrature,
                        vector<double> &weight)
{
    shared_ptr<Basis_Function> basis = basis_function(i);

    
}

void Weight_Function::
calculate_integrals_1d()
{
    for (int i = 0; i < number_of_basis_functions_; ++i)
    {
        vector<double> quadrature;
        vector<double> weight;

        get_basis_quadrature(i,
                             quadrature,
                             weight);
        
        for (int 
    }
}
