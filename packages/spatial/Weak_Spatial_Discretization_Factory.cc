#include "Weak_Spatial_Discretization_Factory.hh"

Weak_Spatial_Discretization_Factory::
Weak_Spatial_Discretization_Factory(shared_ptr<Solid_Geometry> solid_geometry,
                                    vector<shared_ptr<Cartesian_Plane> > const &boundary_surfaces):
    solid_geometry_(solid_geometry),
    boundary_surfaces_(boundary_surfaces_)
{
    int dimension = solid_geometry_->dimension();
    limits_.assign(dimension, vector<double>(2, 0));
    for (shared_ptr<Cartesian_Plane> surface : boundary_surfaces)
    {
        int surface_dimension = surface->surface_dimension();
        double position = surface->position();
        double normal = surface->normal();
        
        if (normal < 0)
        {
            limits_[surface_dimension][0] = position;
        }
        else
        {
            limits_[surface_dimension][1] = position;
        }
    }
}

shared_ptr<Weak_Spatial_Discretization>
square_grid(vector<int> const &number_of_points_1d) const
{
    int dimension = solid_geometry_->dimension();

    // Get points in each dimension
    vector<vector<double> > points_1d(dimension);
    int number_of_points = 1;
    for (int d = 0; d < dimension; ++d)
    {
        points_1d[d].resize(number_of_points_1d[d]);

        double length = limits_[d][1] - limits_[d][0];
        double dl = length / (number_of_points_1d[d] - 1);
        for (int i = 0; i < number_of_points_1d[d]; ++i)
        {
            points_1d[d][i] = limits_[d][0] + i * dl;
        }

        number_of_points *= number_of_points_1d[d];
    }
    
    vector<vector<double> > points(number_of_points);

    switch (dimension)
    {
    case 1:
        for (int i = 0; i < number_of_points_1d[0]; ++i)
        {
            points[i] = {points_1d[0][i]};
        }
        break;
    case 2:
        for (int j = 0; j < number_of_points_1d[1]; ++j)
        {
            for (int i = 0; i < number_of_points_1d[0]; ++i)
            {
                int index = i + number_of_points_1d[0] * j;
                points[index] = {points_1d[0][i], points_1d[1][j]};
            }
        }
        break;
    case 3:
        for (int k = 0; k < number_of_points_1d[2]; ++k)
        {
        for (int j = 0; j < number_of_points_1d[1]; ++j)
        {
            for (int i = 0; i < number_of_points_1d[0]; ++i)
            {
                int index = i + number_of_points_1d[0] * (j + number_of_points_1d[1] * k);
                points[index] = {points_1d[0][i], points_1d[1][j], points_1d[2][k]};
            }
        }
        }
        break;
    }
    
}
