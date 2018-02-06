#include "Heat_Transfer_Factory.hh"

#include "Boundary_Source.hh"
#include "Cartesian_Plane.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Energy_Discretization.hh"
#include "Gauss_Legendre_Quadrature.hh"
#include "Heat_Transfer_Data.hh"
#include "Heat_Transfer_Integration.hh"
#include "Heat_Transfer_Solve.hh"
#include "Integration_Mesh.hh"
#include "Material_Factory.hh"
#include "Region.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Factory.hh"
#include "Weight_Function.hh"

using namespace std;

Heat_Transfer_Factory::
Heat_Transfer_Factory()
{
}

void Heat_Transfer_Factory::
get_solid_1d(double length,
             shared_ptr<Solid_Geometry> &solid,
             vector<shared_ptr<Cartesian_Plane> > &cartesian_surfaces) const
{
    int dimension = 1;
    
    // Get placeholder energy and angular discretizations
    shared_ptr<Energy_Discretization> energy
        = make_shared<Energy_Discretization>(1); // number of groups
    shared_ptr<Angular_Discretization> angular
        = make_shared<Gauss_Legendre_Quadrature>(dimension,
                                                 1, // number of moments
                                                 2); // number of ordinates
    
    // Get surfaces
    cartesian_surfaces.resize(2);
    cartesian_surfaces[0] = make_shared<Cartesian_Plane>(0, // index
                                               dimension,
                                               Surface::Surface_Type::BOUNDARY,
                                               0, // surface dimension
                                               0, // position
                                               -1); // normal
    cartesian_surfaces[1] = make_shared<Cartesian_Plane>(1, // index
                                               dimension,
                                               Surface::Surface_Type::BOUNDARY,
                                               0, // surface dimension
                                               length, // position
                                               1); // normal
    
    // Get boundary source
    Boundary_Source::Dependencies deps;
    deps.angular = Boundary_Source::Dependencies::Angular::ISOTROPIC;
    vector<shared_ptr<Boundary_Source> > boundary_sources(1);
    vector<double> boundary_data(1, 0);
    boundary_sources[0]
        = make_shared<Boundary_Source>(0, // index
                                       deps,
                                       angular,
                                       energy,
                                       boundary_data,
                                       boundary_data); // alpha
    for (shared_ptr<Cartesian_Plane> &surface : cartesian_surfaces)
    {
        surface->set_boundary_source(boundary_sources[0]);
    }
    vector<shared_ptr<Surface> > surfaces;
    for (shared_ptr<Surface> surface : cartesian_surfaces)
    {
        surfaces.push_back(surface);
    }
    // Get material
    Material_Factory mat_factory(angular,
                                 energy); 
    vector<shared_ptr<Material> > materials(1);
    materials[0]
        = mat_factory.get_standard_material(0, // index
                                            {1}, // sigma_t
                                            {0}, // sigma_s
                                            {0}, // nu
                                            {0}, // sigma_f
                                            {0}, // chi
                                            {0}); // source
    
    // Get region
    vector<Surface::Relation> surface_relations(2, Surface::Relation::NEGATIVE);
    vector<shared_ptr<Region> > regions(1);
    regions[0] = make_shared<Region>(0, // index
                                     materials[0],
                                     surface_relations,
                                     surfaces);

    solid =  make_shared<Constructive_Solid_Geometry>(dimension,
                                                      surfaces,
                                                      regions,
                                                      materials,
                                                      boundary_sources);
}

shared_ptr<Weak_Spatial_Discretization> Heat_Transfer_Factory::
get_spatial_discretization_1d(int number_of_points,
                              double radius_num_intervals,
                              double length,
                              bool basis_mls,
                              bool weight_mls,
                              string basis_type,
                              string weight_type) const
{
    int dimension = 1;

    shared_ptr<Solid_Geometry> solid;
    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces;
    get_solid_1d(length,
                 solid,
                 boundary_surfaces);

    shared_ptr<Weight_Function_Options> weight_options
        = make_shared<Weight_Function_Options>();
    weight_options->tau_const = 0.0;
    shared_ptr<Weak_Spatial_Discretization_Options> weak_options
        = make_shared<Weak_Spatial_Discretization_Options>();
    weak_options->integration_ordinates = 2;
    weak_options->include_supg = false;
    
    Weak_Spatial_Discretization_Factory factory(solid,
                                                boundary_surfaces);
    return factory.get_simple_discretization(number_of_points,
                                             radius_num_intervals,
                                             basis_mls,
                                             weight_mls,
                                             basis_type,
                                             weight_type,
                                             weight_options,
                                             weak_options);
}

shared_ptr<Heat_Transfer_Solve> Heat_Transfer_Factory::
get_cylindrical_1d(int number_of_points,
                   double radius_num_intervals,
                   double radius,
                   bool basis_mls,
                   bool weight_mls,
                   string basis_type,
                   string weight_type,
                   std::shared_ptr<Heat_Transfer_Data> data) const
{
    shared_ptr<Weak_Spatial_Discretization> spatial
        = get_spatial_discretization_1d(number_of_points,
                                        radius_num_intervals,
                                        radius,
                                        basis_mls,
                                        weight_mls,
                                        basis_type,
                                        weight_type);
    shared_ptr<Heat_Transfer_Integration_Options> options
        = make_shared<Heat_Transfer_Integration_Options>();
    options->geometry = Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_1D;
    shared_ptr<Integration_Mesh_Options> integration_options
        = make_shared<Integration_Mesh_Options>();
    integration_options->initialize_from_weak_options(spatial->options());
    integration_options->integration_ordinates = 64;
    shared_ptr<Heat_Transfer_Integration> integration
        = make_shared<Heat_Transfer_Integration>(options,
                                                 integration_options,
                                                 data,
                                                 spatial);
    
    return make_shared<Heat_Transfer_Solve>(integration,
                                            spatial);
}
