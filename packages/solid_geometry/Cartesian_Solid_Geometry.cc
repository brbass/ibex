#include "Cartesian_Solid_Geometry.hh"

#include <cmath>

#include "Check.hh"

using namespace std;

Cartesian_Solid_Geometry::
Cartesian_Solid_Geometry(int dimension,
                         vector<int> const &number_of_cells,
                         vector<double> const &total_length,
                         vector<shared_ptr<Material> > const &materials):
    dimension_(dimension),
    number_of_cells_(number_of_cells),
    total_length_(total_length),
    materials_(materials)
{
    cell_length_.assign(dimension, 0);
    for (int d = 0; d < dimension; ++d)
    {
        cell_length_[d] = total_length_[d] / number_of_cells_[d];
    }
}

int Cartesian_Solid_Geometry::
find_region(vector<double> const &position) const
{
    int multiplication_factor = 1;
    int region = 0;
    for (int d = 0; d < dimension_; ++d)
    {
        int c = floor(position[d] / cell_length_[d]);

        if (c < 0 || c > number_of_cells_[d])
        {
            return Geometry_Errors::NO_REGION;
        }
        
        region += c * multiplication_factor;
        multiplication_factor *= number_of_cells_[d];
    }

    return region;
}

int Cartesian_Solid_Geometry::
find_surface(vector<double> const &position) const
{
    AssertMsg(false, "find_surface not yet implemented");
}

int Cartesian_Solid_Geometry::
next_intersection(vector<double> const &initial_position,
                  vector<double> const &initial_direction,
                  int &final_region,
                  double &distance,
                  vector<double> &final_position) const
{
    AssertMsg(false, "next_intersection not yet implemented");
}

int Cartesian_Solid_Geometry::
next_boundary(vector<double> const &initial_position,
              vector<double> const &initial_direction,
              int &boundary_region,
              double &distance,
              vector<double> &final_position) const
{
    AssertMsg(false, "next_boundary not yet implemented");
}

void Cartesian_Solid_Geometry::
optical_distance(vector<double> const &initial_position,
                 vector<double> const &final_position,
                 vector<double> &optical_distance) const
{
    AssertMsg(false, "next_boundary not yet implemented");
}

shared_ptr<Material> Cartesian_Solid_Geometry::
material(vector<double> const &position) const
{
    int region = find_region(position);

    if (region == Geometry_Errors::NO_REGION)
    {
        return shared_ptr<Material>();
    }

    return materials_[region];
}

void Cartesian_Solid_Geometry::
check_class_invariants() const
{
}

void Cartesian_Solid_Geometry::
output(pugi::xml_node &output_node) const
{
}
