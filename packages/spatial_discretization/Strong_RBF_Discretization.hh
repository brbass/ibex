#ifndef Strong_RBF_Discretization_hh
#define Strong_RBF_Discretization_hh

#include "RBF_Point.hh"
#include "Spatial_Discretization.hh"

class Angular_Discretization;
class Energy_Discretization;
class Constructive_Solid_Geometry;
template<class Ordinal, class Scalar> class Symmetric_Sparse_Storage;

class Strong_RBF_Discretization : public Spatial_Discretization
{
public:
    
    Strong_RBF_Discretization(int dimension,
                              int number_of_points,
                              int number_of_internal_points,
                              int number_of_boundary_points,
                              int number_of_neighbors,
                              std::vector<int> const &internal_points,
                              std::vector<int> const &boundary_points,
                              std::vector<std::shared_ptr<RBF_Point> > const &rbf_points,
                              std::shared_ptr<Constructive_Solid_Geometry> solid_geometry);

    // Spatial_Discretization functions
    virtual int dimension() const override
    {
        return dimension_;
    }
    virtual int number_of_points() const override
    {
        return number_of_points_;
    }
    virtual int number_of_internal_points() const
    {
        return number_of_internal_points_;
    }
    virtual int number_of_boundary_points() const override
    {
        return number_of_boundary_points_;
    }
    virtual int number_of_nodes() const override
    {
        return 1;
    }
    virtual std::vector<int> const &boundary_points() const override
    {
        return boundary_points_;
    }
    virtual std::vector<int> const &internal_points() const
    {
        return internal_points_;
    }
    virtual std::shared_ptr<Point> point(int point_index) const override
    {
        return rbf_points_[point_index];
    }
    virtual void output(XML_Node output_node) const override;

    virtual void check_class_invariants() const override;
    
    // Strong_RBF_Discretization functions
    virtual int number_of_neighbors() const
    {
        return number_of_neighbors_;
    }
    virtual std::shared_ptr<RBF_Point> rbf_point(int point_index) const
    {
        return rbf_points_[point_index];
    }
    virtual std::shared_ptr<Constructive_Solid_Geometry> solid_geometry() const
    {
        return solid_geometry_;
    }
    
private:

    int dimension_;
    int number_of_points_;
    int number_of_internal_points_;
    int number_of_boundary_points_;
    int number_of_neighbors_;
    std::vector<int> internal_points_;
    std::vector<int> boundary_points_;
    std::vector<std::shared_ptr<RBF_Point> > rbf_points_;
    std::shared_ptr<Constructive_Solid_Geometry> solid_geometry_;
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;
};

#endif
