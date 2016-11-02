#ifndef Region_hh
#define Region_hh

#include <memory>
#include <vector>

#include "Surface.hh"

class Material;

using std::shared_ptr;
using std::vector;

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

    // Constructor
    Region(int index,
           shared_ptr<Material> material,
           vector<Surface::Relation> const &surface_relations,
           vector<shared_ptr<Surface> > const &surfaces);

    int index() const
    {
        return index_;
    }
    int number_of_surfaces() const
    {
        return surfaces_.size();
    }
    shared_ptr<Material> material() const
    {
        return material_;
    }
    Surface::Relation surface_relation(int s) const
    {
        return surface_relations_[s];
    }
    shared_ptr<Surface> const &surface(int s) const
    {
        return surfaces_[s];
    }
    vector<shared_ptr<Surface> > const &surfaces() const
    {
        return surfaces_;
    }
    
    Relation relation(vector<double> const &point) const;

    virtual void check_class_invariants() const;

    virtual void output(pugi::xml_node &output_node) const;
    
private:

    int index_;
    shared_ptr<Material> material_;
    vector<Surface::Relation> surface_relations_;
    vector<shared_ptr<Surface> > surfaces_;
};

#endif
