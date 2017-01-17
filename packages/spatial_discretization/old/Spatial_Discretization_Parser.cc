#include "Spatial_Discretization_Parser.hh"

#include "Boundary_Source.hh"
#include "Cartesian_Distance.hh"
#include "Cartesian_Overlay.hh"
#include "Derivative_RBF_Function.hh"
#include "Distance.hh"
#include "Energy_Discretization.hh"
#include "Gaussian_RBF.hh"
#include "Inverse_Multiquadric_RBF.hh"
#include "Material.hh"
#include "Multiquadric_RBF.hh"
#include "Optical_Distance.hh"
#include "RBF.hh"
#include "RBF_Discretization.hh"
#include "RBF_Function.hh"
#include "RBF_Point.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Constructive_Solid_Geometry_Parser.hh"
#include "Spatial_Parser_Functions.hh"
#include "Vector_Functions_2D.hh"
#include "Vector_Functions_3D.hh"
#include "Wendland_RBF.hh"

using namespace std;
namespace vf2 = Vector_Functions_2D;
namespace vf3 = Vector_Functions_3D;
namespace spf = Spatial_Parser_Functions;

Spatial_Discretization_Parser::
Spatial_Discretization_Parser(pugi::xml_node &input_file,
                              vector<shared_ptr<Material> > const &materials,
                              vector<shared_ptr<Boundary_Source> > const &boundary_sources):
    Parser(input_file),
    materials_(materials),
    boundary_sources_(boundary_sources)
{
    pugi::xml_node spatial = input_file.child("spatial_discretization");
    
    string spatial_type = XML_Functions::child_value<string>(spatial, "type");
    
    if (spatial_type == "rbf")
    {
        spatial_ = get_rbf(spatial);
    }
    else 
    {
        AssertMsg(false, "spatial discretization type \"" + spatial_type + "\" not found");
    }
}

shared_ptr<RBF_Discretization> Spatial_Discretization_Parser::
get_rbf(pugi::xml_node &spatial)
{
    Constructive_Solid_Geometry_Parser solid_geometry_parser(spatial,
                                                              materials_,
                                                              boundary_sources_);
    
    shared_ptr<Constructive_Solid_Geometry> solid_geometry = solid_geometry_parser.get_ptr();
    
    int dimension = solid_geometry->dimension();
    int number_of_points;
    int number_of_boundary_points;
    int number_of_internal_points;
    vector<int> boundary_points;
    vector<int> internal_points;
    vector<int> surface;
    vector<int> region;
    vector<vector<double> > positions;
    vector<vector<double> > boundary_normal;
    vector<shared_ptr<Material> > material;
    vector<shared_ptr<Boundary_Source> > boundary_source;
    shared_ptr<Cartesian_Overlay> cartesian_overlay;

    int max_attempts = XML_Functions::child_value<int>(spatial, "max_attempts");
    double min_distance_boundary = XML_Functions::child_value<double>(spatial, "min_distance_boundary");
    double min_distance_internal = XML_Functions::child_value<double>(spatial, "min_distance_internal");
    double bounding_radius = XML_Functions::child_value<double>(spatial, "bounding_radius");
    vector<double> bounding_origin = XML_Functions::child_vector<double>(spatial, "bounding_origin", dimension);
    shared_ptr<Distance> distance_metric = make_shared<Cartesian_Distance>(dimension);
    
    spf::get_random_points(solid_geometry,
                           distance_metric,
                           dimension,
                           max_attempts,
                           min_distance_boundary,
                           min_distance_internal,
                           bounding_radius,
                           bounding_origin,
                           number_of_points,
                           number_of_boundary_points,
                           number_of_internal_points,
                           surface,
                           region,
                           boundary_points,
                           internal_points,
                           positions,
                           boundary_normal,
                           material,
                           boundary_source,
                           cartesian_overlay);
    
    string basis_str = XML_Functions::child_value<string>(spatial, "basis_type");
    shared_ptr<RBF> rbf;
    
    if (basis_str == "gaussian")
    {
        rbf = make_shared<Gaussian_RBF>();
    }
    else if (basis_str == "multiquadric")
    {
        rbf = make_shared<Multiquadric_RBF>();
    }
    else if (basis_str == "inverse_multiquadric")
    {
        rbf = make_shared<Inverse_Multiquadric_RBF>();
    }
    else if (basis_str == "wendland30")
    {
        rbf = make_shared<Wendland_RBF>(0);
    }
    else if (basis_str == "wendland31")
    {
        rbf = make_shared<Wendland_RBF>(1);
    }
    else if (basis_str == "wendland32")
    {
        rbf = make_shared<Wendland_RBF>(2);
    }
    else if (basis_str == "wendland33")
    {
        rbf = make_shared<Wendland_RBF>(3);
    }
    else
    {
        AssertMsg(false, "basis_type \"" + basis_str + "\" not found");
    }

    string distance_str = XML_Functions::child_value<string>(spatial, "distance_type");
    shared_ptr<Distance> distance;
    
    if (distance_str == "cartesian")
    {
        distance = make_shared<Cartesian_Distance>(dimension);
    }
    else if (distance_str == "optical")
    {
        distance = make_shared<Optical_Distance>(dimension,
                                                 solid_geometry);
    }
    else
    {
        AssertMsg(false, "distance_type \"" + distance_str + "\" not found");
    }

    string function_str = XML_Functions::child_value<string>(spatial, "function_type");
    shared_ptr<RBF_Function> rbf_function;
    
    if (function_str == "standard")
    {
        rbf_function = make_shared<RBF_Function>(rbf,
                                                 distance);
    }
    else if (function_str == "derivative")
    {
        rbf_function = make_shared<Derivative_RBF_Function>(rbf,
                                                            distance);
    }
    else
    {
        AssertMsg(false, "function_type \"" + function_str + "\" not found");
    }
    
    int number_of_neighbors = XML_Functions::child_value<int>(spatial, "number_of_neighbors");
    if (number_of_neighbors > number_of_points)
    {
        number_of_neighbors = number_of_points;
    }
    
    int number_to_average = XML_Functions::child_value<int>(spatial, "number_to_average");
    double shape_multiplier = XML_Functions::child_value<double>(spatial, "shape_multiplier");
    string shape_method_str
        = XML_Functions::child_value<string>(spatial, "shape_method");

    enum class Shape_Method
    {
        MEAN_DISTANCE,
        MIN_DISTANCE
    };
    
    Shape_Method shape_method;
    if (shape_method_str == "mean_distance")
    {
        shape_method = Shape_Method::MEAN_DISTANCE;
    }
    else if (shape_method_str == "min_distance")
    {
        shape_method = Shape_Method::MIN_DISTANCE;
    }
    
    vector<shared_ptr<RBF_Point> > rbf_points(number_of_points);
    
    int number_of_groups = materials_[0]->energy_discretization()->number_of_groups();
    
    for (int b = 0; b < number_of_boundary_points; ++b)
    {
        int i = boundary_points[b];

        vector<double> position = positions[i];
        vector<double> normal = boundary_normal[b];
        
        vector<int> neighbor_indices;
        vector<double> shape_parameter;
        vector<double> mean_distance;
        
        spf::get_neighbor_information(i,
                                      dimension,
                                      number_of_points,
                                      number_of_groups,
                                      number_of_neighbors,
                                      number_to_average,
                                      distance,
                                      positions,
                                      neighbor_indices,
                                      mean_distance);
        
        switch(shape_method)
        {
        case Shape_Method::MEAN_DISTANCE:
            spf::get_shape_parameter(number_of_groups,
                                     shape_multiplier,
                                     mean_distance,
                                     shape_parameter);
            break;
        case Shape_Method::MIN_DISTANCE:
            if (distance->energy_dependent())
            {
                spf::get_optical_min_distance(number_of_groups,
                                              min_distance_boundary,
                                              material[i]->sigma_t(),
                                              mean_distance);
                
                spf::get_shape_parameter(number_of_groups,
                                         shape_multiplier,
                                         mean_distance,
                                         shape_parameter);
            }
            else
            {
                mean_distance.assign(number_of_groups, min_distance_boundary);
                
                spf::get_shape_parameter(number_of_groups,
                                         shape_multiplier,
                                         mean_distance,
                                         shape_parameter);
            }
        }
        
        rbf_points[i] = make_shared<RBF_Point>(i,
                                               dimension,
                                               number_of_neighbors,
                                               material[i],
                                               boundary_source[b],
                                               rbf_function,
                                               neighbor_indices,
                                               position,
                                               shape_parameter,
                                               mean_distance,
                                               normal);
    }

    for (int p = 0; p < number_of_internal_points; ++p)
    {
        int i = internal_points[p];
        
        vector<double> position = positions[i];
        
        vector<int> neighbor_indices;
        vector<double> shape_parameter;
        vector<double> mean_distance;
        
        spf::get_neighbor_information(i,
                                      dimension,
                                      number_of_points,
                                      number_of_groups,
                                      number_of_neighbors,
                                      number_to_average,
                                      distance,
                                      positions,
                                      neighbor_indices,
                                      mean_distance);
        
        switch(shape_method)
        {
        case Shape_Method::MEAN_DISTANCE:
            spf::get_shape_parameter(number_of_groups,
                                     shape_multiplier,
                                     mean_distance,
                                     shape_parameter);
            break;
        case Shape_Method::MIN_DISTANCE:
            if (distance->energy_dependent())
            {
                spf::get_optical_min_distance(number_of_groups,
                                              min_distance_internal,
                                              material[i]->sigma_t(),
                                              mean_distance);
                
                spf::get_shape_parameter(number_of_groups,
                                         shape_multiplier,
                                         mean_distance,
                                         shape_parameter);
            }
            else
            {
                mean_distance.assign(number_of_groups, min_distance_internal);
                
                spf::get_shape_parameter(number_of_groups,
                                         shape_multiplier,
                                         mean_distance,
                                         shape_parameter);
            }
        }
        
        rbf_points[i] = make_shared<RBF_Point>(i,
                                               dimension,
                                               number_of_neighbors,
                                               material[i],
                                               rbf_function,
                                               neighbor_indices,
                                               position,
                                               shape_parameter,
                                               mean_distance);
    }

    bool store_distances = XML_Functions::child_value<int>(spatial, "store_distances");
    
    return make_shared<RBF_Discretization>(store_distances,
                                           dimension,
                                           number_of_points,
                                           number_of_internal_points,
                                           number_of_boundary_points,
                                           number_of_neighbors,
                                           internal_points,
                                           boundary_points,
                                           rbf_points,
                                           solid_geometry);
}

