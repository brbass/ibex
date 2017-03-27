#include "Weight_Function_Integration.hh"

#include "Basis_Function.hh"
#include "KD_Tree.hh"
#include "Meshless_Function_Factory.hh"
#include "Weight_Function.hh"

Weight_Function_Integration::
Weight_Function_Integration(vector<int> number_of_intervals,
                            shared_ptr<Constructive_Solid_Geometry> solid_geometry,
                            std::vector<std::shared_ptr<Cartesian_Plane> > const &boundaries):
    dimension_(solid_geometry->dimension()),
    number_of_intervals_(number_of_intervals),
    solid_geometry_(solid_geometry),
    boundaries_(boundaries)
{
    Assert(number_of_intervals_.size() == dimension_);
    Assert(solid_geometry_);
    Assert(limits_.size() == dimension_);

    create_mesh();
}

void Weight_Function_Integration::
create_mesh()
{
    switch (dimension)
    {
    case 1:
        return create_mesh_1d();
    case 2:
        return create_mesh_2d();
    default:
        AssertMsg(false, "dimension (" + to_string(dimension) + ") not supported");
    }
}

void create_mesh_1d()
{
    number_of_cells_ = number_of_intervals_[0];
    number_of_nodes_ = number_of_cells + 1;
    number_of_surfaces_ = 2;
    
    cells_.resize(number_of_cells_);
    nodes_.resize(number_of_nodes_);
    surfaces_.resize(number_of_surfaces_);
    
    for (int i = 0; i < number_of_cells_; ++i)
    {
        
    }
}

void create_mesh_2d()
{
    number_of_cells_ = number_of_intervals_[0] * number_of_intervals_[1];
    number_of_nodes_ = (number_of_intervals_[0] + 1) * (number_of_intervals_[1] + 1);
    number_of_surfaces_ = 2 * (number_of_intervals_[0] + number_of_intervals_[1]);
    
    cells_.resize(number_of_cells_);
    nodes_.resize(number_of_nodes_);
    surfaces_.resize(number_of_surfaces_);
    
    
}

void Weight_Function_Integration::
get_weight_functions(int number_of_points,
                     std::vector<Weight_Function::Options> const &options,
                     std::vector<std::vector<int> > const &neighbors,
                     std::vector<std::shared_ptr<Basis_Function> > const &bases,
                     std::vector<std::shared_ptr<Weight_Function> > &weights)
{
}
