#ifndef RBF_Point_hh
#define RBF_Point_hh

#include "Point.hh"
#include "RBF_Function.hh"

// class Epetra_SerialDenseMatrix;
// class Epetra_SerialDenseSolver;

class RBF_Point : public Point
{
public:
    
    RBF_Point(int index,
              int dimension,
              int number_of_neighbors,
              shared_ptr<Material> material,
              shared_ptr<RBF_Function> rbf_function,
              vector<int> const &neighbor_indices,
              vector<double> const &position,
              vector<double> const &shape_parameter,
              vector<double> const &mean_distance);

    RBF_Point(int index,
              int dimension,
              int number_of_neighbors,
              shared_ptr<Material> material,
              shared_ptr<Boundary_Source> boundary_source,
              shared_ptr<RBF_Function> rbf_function,
              vector<int> const &neighbor_indices,
              vector<double> const &position,
              vector<double> const &shape_parameter,
              vector<double> const &mean_distance,
              vector<double> const &normal);
    
    // Point functions
    virtual void output(pugi::xml_node &output_node) const override;
    virtual void check_class_invariants() const override;
    
    // RBF_Point functions
    virtual void set_shape_multiplier(double shape_multiplier);
    virtual void set_neighbors(vector<shared_ptr<RBF_Point> > neighbor_points);
    
    virtual int number_of_neighbors() const
    {
        return number_of_neighbors_;
    }
    virtual vector<int> const &neighbor_indices() const
    {
        return neighbor_indices_;
    }
    virtual vector<double> const &shape_parameter() const
    {
        return shape_parameter_;
    }
    virtual shared_ptr<RBF_Function> rbf_function() const
    {
        return rbf_function_;
    }
    
private:
    
    bool neighbors_set_;
    int number_of_neighbors_;
    vector<int> neighbor_indices_;
    vector<double> shape_parameter_;
    vector<double> mean_distance_;
    vector<shared_ptr<RBF_Point> > neighbor_points_;
    shared_ptr<RBF_Function> rbf_function_;
};

#endif
