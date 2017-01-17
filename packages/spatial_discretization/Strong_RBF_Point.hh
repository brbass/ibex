#ifndef Strong_RBF_Point_hh
#define Strong_RBF_Point_hh

#include "Point.hh"

class RBF_Function;
class XML_Node;

class Strong_RBF_Point : public Point
{
public:
    
    Strong_RBF_Point(int index,
                     int dimension,
                     int number_of_neighbors,
                     std::shared_ptr<Material> material,
                     std::shared_ptr<RBF_Function> rbf_function);

    Strong_RBF_Point(int index,
                     int dimension,
                     int number_of_neighbors,
                     std::shared_ptr<Material> material,
                     std::shared_ptr<Boundary_Source> boundary_source,
                     std::shared_ptr<RBF_Function> rbf_function,
                     std::vector<double> const &normal);
    
    // Point functions
    virtual int index() const override
    {
        return index_;
    }
    virtual int dimension() const override
    {
        return dimension_;
    }
    virtual Point_Type point_type() const override
    {
        return point_type_;
    }
    virtual std::shared_ptr<Material> material() const override
    {
        return material_;
    }
    virtual std::vector<double> const &position() const override
    {
        return position_;
    }
    virtual void output(XML_Node output_node) const override;
    virtual void check_class_invariants() const override;
    
    // Strong_RBF_Point functions
    virtual void set_shape_multiplier(double shape_multiplier);
    virtual void set_neighbors(std::vector<std::shared_ptr<Strong_RBF_Point> > neighbor_points);
    
    virtual int number_of_neighbors() const
    {
        return number_of_neighbors_;
    }
    virtual std::shared_ptr<RBF_Function> rbf_function() const
    {
        return rbf_function_;
    }
    
private:

    bool neighbors_set_;
    int index_;
    int dimension_;
    int number_of_neighbors_;
    Point_Type point_type_;
    std::shared_ptr<Material> material_;
    std::shared_ptr<Boundary_Source> boundary_source_;
    std::shared_ptr<RBF_Function> rbf_function_;
    std::vector<double> position_;
    std::vector<double> normal_;
    std::vector<std::shared_ptr<Strong_RBF_Point> > neighbor_points_;
};

#endif
