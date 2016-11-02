#include "Region.hh"

#include "Check.hh"
#include "Material.hh"

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
output(pugi::xml_node &output_node) const
{
    // General information
    
    pugi::xml_node node = output_node.append_child("region");
    XML_Functions::append_child(node, index_, "index");
    XML_Functions::append_child(node, material_->index(), "material_index");

    // Surfaces

    pugi::xml_node surfaces_node = node.append_child("surface_relations");
    
    int number_of_surface_relations = surface_relations_.size();

    for (int i = 0; i < number_of_surface_relations; ++i)
    {
        pugi::xml_node surface_node = surfaces_node.append_child("surface_relation");

        XML_Functions::append_child(surface_node, surfaces_[i]->index(), "surface");
        
        switch(surface_relations_[i])
        {
        case Surface::Relation::POSITIVE:
            XML_Functions::append_child(surface_node, "positive", "relation");
            break;
        case Surface::Relation::NEGATIVE:
            XML_Functions::append_child(surface_node, "negative", "relation");
            break;
        case Surface::Relation::EQUAL:
            XML_Functions::append_child(surface_node, "equal", "relation");
            break;
        }
    }
}
