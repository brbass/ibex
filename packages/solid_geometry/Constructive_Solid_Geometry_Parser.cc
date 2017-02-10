#include "Constructive_Solid_Geometry_Parser.hh"

#include <iterator>

#include "Boundary_Source.hh"
#include "Cartesian_Plane.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Cylinder_2D.hh"
#include "Cylinder_3D.hh"
#include "Material.hh"
#include "Plane_1D.hh"
#include "Plane_2D.hh"
#include "Plane_3D.hh"
#include "Region.hh"
#include "Sphere_3D.hh"
#include "Surface.hh"
#include "XML_Node.hh"

using namespace std;

Constructive_Solid_Geometry_Parser::
Constructive_Solid_Geometry_Parser(vector<shared_ptr<Material> > materials,
                                   vector<shared_ptr<Boundary_Source> > boundary_sources):
    materials_(materials),
    boundary_sources_(boundary_sources)
{
}

shared_ptr<Constructive_Solid_Geometry> Constructive_Solid_Geometry_Parser::
parse_from_xml(XML_Node solid_node)
{
    int dimension = solid_node.get_child_value<int>("dimension");
    
    vector<shared_ptr<Surface> > surfaces = get_surfaces(solid_node.get_child("surfaces"),
                                                         dimension);
    
    vector<shared_ptr<Region> > regions = get_regions(solid_node.get_child("regions"),
                                                      surfaces);
    
    return make_shared<Constructive_Solid_Geometry>(dimension,
                                                    surfaces,
                                                    regions,
                                                    materials_,
                                                    boundary_sources_);
}

vector<shared_ptr<Surface> > Constructive_Solid_Geometry_Parser::
get_surfaces(XML_Node surfaces_node,
             int dimension)
{
    int number_of_surfaces = surfaces_node.get_child_value<int>("number_of_surfaces");
    
    vector<shared_ptr<Surface> > surfaces(number_of_surfaces);

    int checksum = 0;
    for (XML_Node surface_node = surfaces_node.get_child("surface");
         surface_node;
         surface_node = surface_node.get_sibling("surface",
                                                 false))
    {
        int index = surface_node.get_attribute<int>("index");
        string shape = surface_node.get_attribute<string>("shape");
        
        shared_ptr<Surface> surface;
        
        Surface::Surface_Type surface_type;
        
        string type = surface_node.get_attribute<string>("type");
        
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

        if (shape == "cartesian_plane")
        {
            int surface_dimension = surface_node.get_child_value<int>("surface_dimension");
            double position = surface_node.get_child_value<double>("position");
            double normal = surface_node.get_child_value<double>("normal");

            surface = make_shared<Cartesian_Plane>(index,
                                                   dimension,
                                                   surface_type,
                                                   surface_dimension,
                                                   position,
                                                   normal);
        }
        else if (shape == "plane")
        {
            vector<double> origin = surface_node.get_child_vector<double>("origin", dimension);
            vector<double> normal = surface_node.get_child_vector<double>("normal", dimension);

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
            double radius = surface_node.get_child_value<double>("radius");
            vector<double> origin = surface_node.get_child_vector<double>("origin", dimension);
            
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
                vector<double> direction = surface_node.get_child_vector<double>("direction", dimension);
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
            double radius = surface_node.get_child_value<double>("radius");
            vector<double> origin = surface_node.get_child_vector<double>("origin", dimension);
            
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
            AssertMsg(false, "surface shape (" + shape + ") not found");
        }
        
        if (surface_type == Surface::Surface_Type::BOUNDARY)
        {
            int boundary_source_number = surface_node.get_child_value<int>("boundary_source");
            
            surface->set_boundary_source(boundary_sources_[boundary_source_number]);
        }
        
        surfaces[index] = surface;

        checksum += index;
    } // surfaces

    int checksum_expected = number_of_surfaces * (number_of_surfaces - 1) / 2;
    AssertMsg(checksum == checksum_expected, "Surface indexing incorrect");
    
    return surfaces;
}

vector<shared_ptr<Region> > Constructive_Solid_Geometry_Parser::
get_regions(XML_Node regions_node,
            vector<shared_ptr<Surface> > const &surfaces)
{
    int number_of_regions = regions_node.get_child_value<int>("number_of_regions");
    
    vector<shared_ptr<Region> > regions(number_of_regions);

    int checksum = 0;
    for (XML_Node region_node = regions_node.get_child("region");
         region_node;
         region_node = region_node.get_sibling("region",
                                               false))
    {
        int index = region_node.get_attribute<int>("index");
        int material_number = region_node.get_attribute<int>("material");
        
        vector<Surface::Relation> surface_relations;
        vector<shared_ptr<Surface> > local_surfaces;
        
        for (XML_Node surface_relation_node = region_node.get_child("surface_relation");
             surface_relation_node;
             surface_relation_node = surface_relation_node.get_sibling("surface_relation",
                                                                       false))
        {
            int relation_number = surface_relation_node.get_attribute<int>("surface");
            string relation_string = surface_relation_node.get_attribute<string>("relation");
            
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

        checksum += index;
    }
    int checksum_expected = number_of_regions * (number_of_regions - 1) / 2;
    AssertMsg(checksum == checksum_expected, "Region indexing incorrect");
    
    return regions;
}

