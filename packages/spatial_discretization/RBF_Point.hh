#ifndef RBF_Point_hh
#define RBF_Point_hh

#include "Point.hh"
#include "RBF_Function.hh"

class XML_Node;

class RBF_Point : public Point
{
public:
    
    RBF_Point(int index,
              int dimension,
              int number_of_neighbors,
              std::shared_ptr<Material> material,
              std::shared_ptr<RBF_Function> rbf_function,
              std::vector<int> const &neighbor_indices,
              std::vector<double> const &position,
              std::vector<double> const &shape_parameter,
              std::vector<double> const &mean_distance);

    RBF_Point(int index,
              int dimension,
              int number_of_neighbors,
              std::shared_ptr<Material> material,
              std::shared_ptr<Boundary_Source> boundary_source,
              std::shared_ptr<RBF_Function> rbf_function,
              std::vector<int> const &neighbor_indices,
              std::vector<double> const &position,
              std::vector<double> const &shape_parameter,
              std::vector<double> const &mean_distance,
              std::vector<double> const &normal);
    
    // Point functions
    virtual void output(XML_Node output_node) const override;
    virtual void check_class_invariants() const override;
    
    // RBF_Point functions
    virtual void set_shape_multiplier(double shape_multiplier);
    virtual void set_neighbors(std::vector<std::shared_ptr<RBF_Point> > neighbor_points);
    
    virtual int number_of_neighbors() const
    {
        return number_of_neighbors_;
    }
    virtual std::vector<int> const &neighbor_indices() const
    {
        return neighbor_indices_;
    }
    virtual std::vector<double> const &shape_parameter() const
    {
        return shape_parameter_;
    }
    virtual std::shared_ptr<RBF_Function> rbf_function() const
    {
        return rbf_function_;
    }
    
private:
    
    bool neighbors_set_;
    int number_of_neighbors_;
    std::vector<int> neighbor_indices_;
    std::vector<double> shape_parameter_;
    std::vector<double> mean_distance_;
    std::vector<std::shared_ptr<RBF_Point> > neighbor_points_;
    std::shared_ptr<RBF_Function> rbf_function_;
};

#endif
