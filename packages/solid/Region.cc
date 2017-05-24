#include "Region.hh"

#include "Check.hh"
#include "Conversion.hh"
#include "Material.hh"
#include "XML_Node.hh"

using namespace std;

Region::
Region(int index,
       shared_ptr<Material> material,
       vector<Surface::Relation> const &surface_relations,
       vector<shared_ptr<Surface> > const &surfaces):
    index_(index),
    material_(material),
    surface_relations_(surface_relations),
    surfaces_(surfaces)
{
}

Region::Relation Region::
relation(vector<double> const &point) const
{
    Surface::Relation surface_relation;
    
    for (int i = 0; i < surfaces_.size(); ++i)
    {
        surface_relation = surfaces_[i]->relation(point,
                                                  false);
        
        if (surface_relation != surface_relations_[i])
        {
            return Region::Relation::OUTSIDE;
        }
    }
    
    return Region::Relation::INSIDE;
}

void Region::
check_class_invariants() const
{
    Assert(surface_relations_.size() == surfaces_.size());

    int number_of_surfaces = surfaces_.size();
    
    for (int i = 0; i < number_of_surfaces; ++i)
    {
        Assert(surfaces_[i]);
    }
}

void Region::
output(XML_Node output_node) const
{
    // General information
    
    output_node.set_attribute(index_, "index");
    output_node.set_child_value(material_->index(), "material_index");

    // Surfaces

    XML_Node surfaces_node = output_node.append_child("surface_relations");
    
    int number_of_surface_relations = surface_relations_.size();
    
    for (int i = 0; i < number_of_surface_relations; ++i)
    {
        XML_Node surface_node = surfaces_node.append_child("surface_relation");
        
        surface_node.set_attribute(surfaces_[i]->index(), "surface");
        
        switch(surface_relations_[i])
        {
        case Surface::Relation::POSITIVE:
            surface_node.set_attribute("positive", "relation");
            break;
        case Surface::Relation::NEGATIVE:
            surface_node.set_attribute("negative", "relation");
            break;
        case Surface::Relation::EQUAL:
            surface_node.set_attribute("equal", "relation");
            break;
        }
    }
}

shared_ptr<Conversion<Region::Relation, string> > Region::
relation_conversion() const
{
    vector<pair<Relation, string> > conversions
        = {{Relation::INSIDE, "inside"},
           {Relation::OUTSIDE, "outside"}};
    return make_shared<Conversion<Relation, string> >(conversions);
}
