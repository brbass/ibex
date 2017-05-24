#ifndef Basis_Function_hh 
#define Basis_Function_hh

#include <memory>
#include <vector>

class Cartesian_Plane;
template<class T1, class T2> class Conversion;
class Meshless_Function;
class XML_Node;

class Basis_Function
{
public:

    enum class Point_Type
    {
        INTERNAL,
        BOUNDARY
    };
    std::shared_ptr<Conversion<Point_Type, std::string> > point_type_conversion() const;
    
    // Constructor
    Basis_Function(int index,
                   int dimension,
                   std::shared_ptr<Meshless_Function> meshless_function,
                   std::vector<std::shared_ptr<Cartesian_Plane> > boundary_surfaces);
    
    // Data access
    virtual int index() const
    {
        return index_;
    }
    virtual int boundary_index() const
    {
        return boundary_index_;
    }
    virtual int dimension() const
    {
        return dimension_;
    }
    virtual int number_of_boundary_surfaces() const
    {
        return number_of_boundary_surfaces_;
    }
    virtual double radius() const
    {
        return radius_;
    }
    virtual Point_Type point_type() const
    {
        return point_type_;
    }
    virtual std::vector<double> const &position() const
    {
        return position_;
    }
    virtual std::shared_ptr<Meshless_Function> function() const
    {
        return meshless_function_;
    }
    virtual std::shared_ptr<Cartesian_Plane> boundary_surface(int i) const
    {
        return boundary_surfaces_[i];
    }
    virtual void output(XML_Node output_node) const;
    virtual void check_class_invariants() const;
    
    // Set data 
    virtual void set_boundary_index(int index)
    {
        boundary_index_ = index;
    }
    
private:
    
    // Data
    int index_;
    int boundary_index_;
    int dimension_;
    int number_of_boundary_surfaces_;
    double radius_;
    Point_Type point_type_;
    std::vector<double> position_;
    std::shared_ptr<Meshless_Function> meshless_function_;
    std::vector<std::shared_ptr<Cartesian_Plane> > boundary_surfaces_;
};

#endif
