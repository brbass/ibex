#ifndef Region_hh
#define Region_hh

#include <memory>
#include <vector>

#include "Surface.hh"

template<class T1, class T2> class Conversion;
class Material;
class XML_Node;

/* 
   Describes a solid geometry region
   
   The memory for the regions must be allocated before the data is initialized,
   as a region can be defined in terms of other regions. Because of this, the
   "initialize" function must be called after the constructor. This allows 
   construction of the regions in any order without regard for interdependencies.
*/
class Region
{
public:

    // Is the point inside or outside of the region
    enum class Relation
    {
        INSIDE,
        OUTSIDE
    };

    std::shared_ptr<Conversion<Relation, std::string> > relation_conversion() const;
    
    // Constructor
    Region(int index,
           std::shared_ptr<Material> material,
           std::vector<Surface::Relation> const &surface_relations,
           std::vector<std::shared_ptr<Surface> > const &surfaces);

    int index() const
    {
        return index_;
    }
    int number_of_surfaces() const
    {
        return surfaces_.size();
    }
    std::shared_ptr<Material> material() const
    {
        return material_;
    }
    Surface::Relation surface_relation(int s) const
    {
        return surface_relations_[s];
    }
    std::shared_ptr<Surface> const &surface(int s) const
    {
        return surfaces_[s];
    }
    std::vector<std::shared_ptr<Surface> > const &surfaces() const
    {
        return surfaces_;
    }
    
    Relation relation(std::vector<double> const &point) const;

    virtual void check_class_invariants() const;

    virtual void output(XML_Node output_node) const;
    
private:
    
    int index_;
    std::shared_ptr<Material> material_;
    std::vector<Surface::Relation> surface_relations_;
    std::vector<std::shared_ptr<Surface> > surfaces_;
};

#endif
