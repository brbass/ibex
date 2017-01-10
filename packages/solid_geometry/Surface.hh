#ifndef Surface_hh
#define Surface_hh

#include <memory>
#include <vector>

#include "Check.hh"

class Boundary_Source;
class XML_Node;

/*
  General class for a solid geometry surface
*/
class Surface
{
public:

    /* Class of surface */
    enum class Surface_Class
    {
        CARTESIAN_PLANE,
        PLANE,
        SPHERE,
        CYLINDER    
    };

    /* Relationship of point to surface */
    enum class Relation
    {
        POSITIVE,
        NEGATIVE,
        EQUAL,
        OUTSIDE = POSITIVE,
        INSIDE = NEGATIVE
    };
    
    /* Classification of surface */
    enum class Surface_Type
    {
        BOUNDARY,
        INTERNAL
    };
    
    /* Types of intersection of vector and surface */
    enum class Intersection
    {
        INTERSECTS, // has intersection
        PARALLEL, // has no intersection, parallel to surface
        NONE, // has no intersection, negative or positive
        NEGATIVE, // only has negative intersection
        TANGEANT // only intersects at one infinitesimal point
    };

    /* Constructor */
    Surface(int index,
            int dimension,
            Surface_Type surface_type);
    
    /* Surface index */
    virtual int index()
    {
        return index_;
    }
    
    /* Number of spatial dimensions */
    virtual int dimension()
    {
        return dimension_;
    }
    
    /* Returns base class of surface */
    virtual Surface_Class surface_class() const = 0;
    
    /* Returns surface type */
    virtual Surface_Type surface_type() const
    {
        return surface_type_;
    }

    /* Returns relationship between point and surface */
    virtual Relation relation(std::vector<double> const &position,
                              bool check_equality = false) const = 0;
    
    /* Returns distance to the surface */
    virtual double distance(std::vector<double> const &position) const
    {
        AssertMsg(false, "Not yet implemented for this surface");
        
        return 0;
    }

    /* Type of intersection of streaming particle with surface
       If type is TANGEANT or PARALLEL, the distance and position are
       returned. Otherwise, the distance and position remain unchanged.*/
    virtual Intersection intersection(std::vector<double> const &initial_position,
                                      std::vector<double> const &initial_direction,
                                      double &distance,
                                      std::vector<double> &final_position) const = 0;

    /* Normal direction of surface at a point on the surface */
    virtual bool normal_direction(std::vector<double> const &position,
                                  std::vector<double> &normal,
                                  bool check_normal = true) const = 0;

    /* Reflected direction */
    virtual bool reflected_direction(std::vector<double> const &position,
                                     std::vector<double> const &initial_direction,
                                     std::vector<double> &final_direction,
                                     bool check_normal = true);

    /* Boundary source (for boundary surfaces) */
    virtual std::shared_ptr<Boundary_Source> boundary_source()
    {
        return boundary_source_;
    }
    
    virtual void set_boundary_source(std::shared_ptr<Boundary_Source> boundary_source)
    {
        boundary_source_ = boundary_source;
    }

    virtual void check_class_invariants() const = 0;
    
    virtual void output(XML_Node output_node) const = 0;
    
protected:

    int index_;
    int dimension_;
    Surface_Type surface_type_;
    double relation_tolerance_;
    double intersection_tolerance_;
    double normal_tolerance_;
    
    std::shared_ptr<Boundary_Source> boundary_source_;
};

#endif
