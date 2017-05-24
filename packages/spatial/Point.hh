#ifndef Point_hh
#define Point_hh

#include <memory>
#include <string>
#include <vector>

class Boundary_Source;
template<class T1, class T2> class Conversion;
class Material;
class XML_Node;

class Point
{
public:
    
    enum class Point_Type
    {
        INTERNAL,
        BOUNDARY
    };
    std::shared_ptr<Conversion<Point_Type, std::string> > point_type_conversion() const;

    // Constructor
    Point();

    // Data access
    virtual int index() const = 0;
    virtual int dimension() const = 0;
    virtual int number_of_nodes() const = 0;
    virtual Point_Type point_type() const = 0;
    virtual std::shared_ptr<Material> material() const = 0;
    virtual std::vector<double> const &position() const = 0;
    
    // Data output and checking
    virtual void output(XML_Node output_node) const = 0;
    virtual void check_class_invariants() const = 0;
};

#endif
