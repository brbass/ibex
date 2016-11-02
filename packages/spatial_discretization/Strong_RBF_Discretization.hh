#ifndef RBF_Discretization_hh
#define RBF_Discretization_hh

#include "RBF_Point.hh"
#include "Spatial_Discretization.hh"

class Angular_Discretization;
class Energy_Discretization;
class Constructive_Solid_Geometry;
template<class Ordinal, class Scalar> class Symmetric_Sparse_Storage;

class RBF_Discretization : public Spatial_Discretization
{
public:
    
    RBF_Discretization(bool store_distances,
                       int dimension,
                       int number_of_points,
                       int number_of_internal_points,
                       int number_of_boundary_points,
                       int number_of_neighbors,
                       vector<int> const &internal_points,
                       vector<int> const &boundary_points,
                       vector<shared_ptr<RBF_Point> > const &rbf_points,
                       shared_ptr<Constructive_Solid_Geometry> solid_geometry);

    // Spatial_Discretization functions
    virtual int dimension() const override
    {
        return dimension_;
    }
    virtual int number_of_points() const override
    {
        return number_of_points_;
    }
    virtual int number_of_internal_points() const override
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
    virtual vector<int> const &boundary_points() const override
    {
        return boundary_points_;
    }
    virtual vector<int> const &internal_points() const override
    {
        return internal_points_;
    }
    virtual shared_ptr<Point> point(int point_index) const override
    {
        return rbf_points_[point_index];
    }
    virtual void output(pugi::xml_node &output_node) const override;

    virtual void check_class_invariants() const override;
    
    // RBF_Discretization functions
    virtual int number_of_neighbors() const
    {
        return number_of_neighbors_;
    }
    virtual shared_ptr<RBF_Point> rbf_point(int point_index) const
    {
        return rbf_points_[point_index];
    }
    virtual shared_ptr<Constructive_Solid_Geometry> solid_geometry() const
    {
        return solid_geometry_;
    }
    virtual shared_ptr<Symmetric_Sparse_Storage<int, double> > distances(int g) const
    {
        return distances_[g];
    }
    virtual double get_distance(int i,
                                int i0,
                                int g) const;
    virtual double basis(int i,
                         int i0,
                         int g,
                         int o) const;
    virtual vector<double> gradient_basis(int i,
                                          int i0,
                                          int g,
                                          int o) const;
    
    virtual void calculate_distances() const;
    
private:

    bool store_distances_;
    int dimension_;
    int number_of_points_;
    int number_of_internal_points_;
    int number_of_boundary_points_;
    int number_of_neighbors_;
    vector<int> internal_points_;
    vector<int> boundary_points_;
    vector<shared_ptr<RBF_Point> > rbf_points_;
    shared_ptr<Constructive_Solid_Geometry> solid_geometry_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;
    mutable vector<shared_ptr<Symmetric_Sparse_Storage<int, double> > > distances_;
};

#endif
