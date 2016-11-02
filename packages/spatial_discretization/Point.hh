#ifndef Point_hh
#define Point_hh

#include <string>
#include <vector>

#include "Boundary_Source.hh"
#include "Material.hh"

using std::string;
using std::vector;

class Point
{
public:

    enum class Point_Type
    {
        INTERNAL,
        BOUNDARY
    };
    
    Point(int index,
          int dimension,
          shared_ptr<Material> material,
          vector<double> const &position);
    
    Point(int index,
          int dimension,
          shared_ptr<Material> material,
          shared_ptr<Boundary_Source> boundary_source,
          vector<double> const &position,
          vector<double> const &normal);
    
    virtual int index() const
    {
        return index_;
    }
    virtual int dimension() const
    {
        return dimension_;
    }
    virtual Point_Type point_type() const
    {
        return point_type_;
    }
    virtual shared_ptr<Material> material() const
    {
        return material_;
    }
    virtual shared_ptr<Boundary_Source> boundary_source() const
    {
        return boundary_source_;
    }
    virtual vector<double> const &position() const
    {
        return position_;
    }
    virtual vector<double> const &normal() const
    {
        return normal_;
    }
    
    virtual void check_class_invariants() const;
    
    virtual void output(pugi::xml_node &output_node) const;
    
protected:

    virtual string point_type_string() const;
    
    int index_;
    int dimension_;
    Point_Type point_type_;
    shared_ptr<Material> material_;
    shared_ptr<Boundary_Source> boundary_source_;
    vector<double> position_;
    vector<double> normal_;
};

#endif
