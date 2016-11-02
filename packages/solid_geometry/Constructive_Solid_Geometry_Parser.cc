#include "Constructive_Solid_Geometry_Parser.hh"

#include <iterator>

#include "Boundary_Source.hh"
#include "Cylinder_2D.hh"
#include "Cylinder_3D.hh"
#include "Material.hh"
#include "Plane_1D.hh"
#include "Plane_2D.hh"
#include "Plane_3D.hh"
#include "Sphere_3D.hh"

using namespace std;

Constructive_Solid_Geometry_Parser::
Constructive_Solid_Geometry_Parser(pugi::xml_node &input_file,
                      vector<shared_ptr<Material> > materials,
                      vector<shared_ptr<Boundary_Source> > boundary_sources):
    Parser(input_file),
    materials_(materials),
    boundary_sources_(boundary_sources)
{
    pugi::xml_node solid_node = input_file.child("solid_geometry");

    solid_ = get_solid(solid_node);
}

shared_ptr<Constructive_Solid_Geometry> Constructive_Solid_Geometry_Parser::
get_solid(pugi::xml_node &solid_node)
{
    int dimension = XML_Functions::child_value<int>(solid_node, "dimension");
    
    pugi::xml_node surfaces_node = solid_node.child("surfaces");
    vector<shared_ptr<Surface> > surfaces = get_surfaces(surfaces_node,
                                                         dimension);

    pugi::xml_node regions_node = solid_node.child("regions");
    vector<shared_ptr<Region> > regions = get_regions(regions_node,
                                                      surfaces);
    
    return make_shared<Constructive_Solid_Geometry>(dimension,
                                                     surfaces,
                                                     regions,
                                                     materials_,
                                                     boundary_sources_);
}

vector<shared_ptr<Surface> > Constructive_Solid_Geometry_Parser::
get_surfaces(pugi::xml_node &surfaces_node,
             int dimension)
{
    int number_of_surfaces = distance(surfaces_node.children("surface").begin(),
                                      surfaces_node.children("surface").end());
    
    vector<shared_ptr<Surface> > surfaces(number_of_surfaces);

    for (pugi::xml_node surface_node = surfaces_node.child("surface"); surface_node; surface_node = surface_node.next_sibling("surface"))
    {
        int index = XML_Functions::child_value<int>(surface_node, "index");
        string shape = XML_Functions::child_value<string>(surface_node, "shape");
        
        shared_ptr<Surface> surface;
        
        Surface::Surface_Type surface_type;
        
        string type = XML_Functions::child_value<string>(surface_node, "type");

        shared_ptr<Boundary_Source> boundary_source;
        
        if (type == "boundary")
        {
            surface_type = Surface::Surface_Type::BOUNDARY;
        }
        else if (type == "internal")
        {
            surface_type = Surface::Surface_Type::INTERNAL;
        }
        else
        {
            AssertMsg(false, "surface type not found");
        }

        if (shape == "plane")
        {
            vector<double> origin = XML_Functions::child_vector<double>(surface_node, "origin", dimension);
            vector<double> normal = XML_Functions::child_vector<double>(surface_node, "normal", dimension);

            switch(dimension)
            {
            case 1:
                surface = make_shared<Plane_1D>(index,
                                                surface_type,
                                                origin,
                                                normal);
                break;
            case 2:
                surface = make_shared<Plane_2D>(index,
                                                surface_type,
                                                origin,
                                                normal);
                break;
            case 3:
                surface = make_shared<Plane_3D>(index,
                                                surface_type,
                                                origin,
                                                normal);
                break;
            default:
                AssertMsg(false, "dimension not found");
            }
        }
        else if (shape == "cylinder")
        {
            double radius = XML_Functions::child_value<double>(surface_node, "radius");
            vector<double> origin = XML_Functions::child_vector<double>(surface_node, "origin", dimension);

            switch(dimension)
            {
            case 2:
                surface = make_shared<Cylinder_2D>(index,
                                                   surface_type,
                                                   radius,
                                                   origin);
                break;
            case 3:
            {
                vector<double> direction = XML_Functions::child_vector<double>(surface_node, "direction", dimension);
                surface = make_shared<Cylinder_3D>(index,
                                                   surface_type,
                                                   radius,
                                                   origin,
                                                   direction);
                break;
            }
            default:                
                AssertMsg(false, "dimension not found");
                
            }
        }
        else if (shape == "sphere")
        {
            double radius = XML_Functions::child_value<double>(surface_node, "radius");
            vector<double> origin = XML_Functions::child_vector<double>(surface_node, "origin", dimension);

            switch(dimension)
            {
            case 3:
                surface = make_shared<Sphere_3D>(index,
                                                 surface_type,
                                                 radius,
                                                 origin);
                break;
            default:
                AssertMsg(false, "dimension not found");
            }
        }
        else
        {
            AssertMsg(false, "surface shape not found");
        }

        if (surface->surface_type() == Surface::Surface_Type::BOUNDARY)
        {
            int boundary_source_number = XML_Functions::child_value<int>(surface_node, "boundary_source");
            
            surface->set_boundary_source(boundary_sources_[boundary_source_number]);
        }
        
        surfaces[index] = surface;
    }
    
    return surfaces;
}

vector<shared_ptr<Region> > Constructive_Solid_Geometry_Parser::
get_regions(pugi::xml_node &regions_node,
            vector<shared_ptr<Surface> > const &surfaces)
{
    int number_of_regions = distance(regions_node.children("region").begin(),
                                     regions_node.children("region").end());
    
    vector<shared_ptr<Region> > regions(number_of_regions);

    for (pugi::xml_node region_node = regions_node.child("region"); region_node; region_node = region_node.next_sibling("region"))
    {
        int index = XML_Functions::child_value<int>(region_node, "index");
        int material_number = XML_Functions::child_value<int>(region_node, "material");
        
        vector<Surface::Relation> surface_relations;
        vector<shared_ptr<Surface> > local_surfaces;
 
        for (pugi::xml_node surface_relation_node = region_node.child("surface_relation"); surface_relation_node; surface_relation_node = surface_relation_node.next_sibling("surface_relation"))
        {
            int relation_number = XML_Functions::child_value<int>(surface_relation_node, "surface");
            string relation_string = XML_Functions::child_value<string>(surface_relation_node, "relation");
            
            Surface::Relation relation;

            if (relation_string == "positive")
            {
                relation = Surface::Relation::POSITIVE;
            }
            else if (relation_string == "negative")
            {
                relation = Surface::Relation::NEGATIVE;
            }
            else if (relation_string == "equal")
            {
                relation = Surface::Relation::EQUAL;
            }
            else if (relation_string == "outside")
            {
                relation = Surface::Relation::OUTSIDE;
            }
            else if (relation_string == "inside")
            {
                relation = Surface::Relation::INSIDE;
            }
            else
            {
                AssertMsg(false, "surface relation not found");
            }

            surface_relations.push_back(relation);
            local_surfaces.push_back(surfaces[relation_number]);
        }       
        
        regions[index] = make_shared<Region>(index,
                                             materials_[material_number],
                                             surface_relations,
                                             local_surfaces);
    }
    
    return regions;
}
