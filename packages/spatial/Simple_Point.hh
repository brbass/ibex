#ifndef Simple_Point_hh
#define Simple_Point_hh

#include <memory>
#include <vector>

#include "Point.hh"

class Material;
class XML_Node;

class Simple_Point : public Point
{
public:

    Simple_Point(int index,
                 int dimension,
                 Point_Type point_type,
                 std::shared_ptr<Material> material,
                 std::vector<double> const &position);

    virtual int index() const override
    {
        return index_;
    }
    virtual int dimension() const override
    {
        return dimension_;
    }
    virtual int number_of_nodes() const override
    {
        return 1;
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
    
    // Data output and checking
    virtual void output(XML_Node output_node) const override;
    virtual void check_class_invariants() const override;

private:

    int index_;
    int dimension_;
    Point_Type point_type_;
    std::shared_ptr<Material> material_;
    std::vector<double> position_;
};

#endif
