#include <iostream>
#include <string>
#include <vector>

#include "Angular_Discretization_Factory.hh"
#include "Boundary_Source.hh"
#include "Cartesian_Plane.hh"
#include "Combined_SUPG_Fission.hh"
#include "Combined_SUPG_Scattering.hh"
#include "Constructive_Solid_Geometry.hh"
#include "Cylinder_2D.hh"
#include "Discrete_Normalization_Operator.hh"
#include "Energy_Discretization.hh"
#include "Fission.hh"
#include "LDFE_Quadrature.hh"
#include "Material.hh"
#include "Material_Factory.hh"
#include "Moment_To_Discrete.hh"
#include "Moment_Weighting_Operator.hh"
#include "Random_Number_Generator.hh"
#include "Region.hh"
#include "SUPG_Fission.hh"
#include "SUPG_Moment_To_Discrete.hh"
#include "SUPG_Scattering.hh"
#include "Scattering.hh"
#include "Vector_Operator_Functions.hh"
#include "Weak_Spatial_Discretization.hh"
#include "Weak_Spatial_Discretization_Factory.hh"

using namespace std;

// Create a three-group, anisotropically-scattering pincell problem with reflective boundaries
// Test whether the various scattering and fission operators return identical results when they should

vector<shared_ptr<Material> >
get_materials(shared_ptr<Angular_Discretization> angular,
              shared_ptr<Energy_Discretization> energy)
{
    int number_of_materials = 2;
    Material_Factory material_factory(angular,
                                      energy);
    vector<shared_ptr<Material> > materials(number_of_materials);
    {
        int index = 0;
        vector<double> sigma_t = {1.0, 2.0, 3.0};
        vector<double> sigma_s = {0.02, 0.05, 0.01,
                                  0.06, 0.58, 0.16,
                                  0.10, 0.69, 0.22,
                                  0.001, 0.02, 0.0001,
                                  0.005, 0.1, 0.005,
                                  0.008, 0.2, 0.02};
        vector<double> nu = {2.4, 1.9, 2.1};
        vector<double> sigma_f = {0.1, 0.2, 0.5};
        vector<double> chi = {0.5, 0.2, 0.3};
        vector<double> internal_source = {0.0, 0.0, 0.0};
        
        materials[0]
            = material_factory.get_standard_material(index,
                                                     sigma_t,
                                                     sigma_s,
                                                     nu,
                                                     sigma_f,
                                                     chi,
                                                     internal_source);
    }
    {
        int index = 0;
        vector<double> sigma_t = {1.1, 2.1, 3.1};
        vector<double> sigma_s = {0.02, 0.1, 0.12,
                                  0.6, 0.8, 0.65,
                                  0.10, 0.09, 0.24,
                                  0.01, 0.002, 0.01,
                                  0.05, 0.012, 0.05,
                                  0.0084, 0.1, 0.07};
        vector<double> nu = {0.0, 0.0, 0.0};
        vector<double> sigma_f = {0.0, 0.0, 0.0};
        vector<double> chi = {0.0, 0.0, 0.0};
        vector<double> internal_source = {0.0, 0.0, 0.0};
        
        materials[1]
            = material_factory.get_standard_material(index,
                                                     sigma_t,
                                                     sigma_s,
                                                     nu,
                                                     sigma_f,
                                                     chi,
                                                     internal_source);
    }
    
    return materials;
 }

vector<shared_ptr<Boundary_Source> >
get_boundary_sources(shared_ptr<Angular_Discretization> angular,
                     shared_ptr<Energy_Discretization> energy)
{
    int number_of_boundary_sources = 1;
    vector<shared_ptr<Boundary_Source> > boundary_sources(number_of_boundary_sources);

    int index = 0;
    Boundary_Source::Dependencies dependencies;
    vector<double> boundary_source = {0.0, 0.0, 0.0};
    vector<double> alpha = {1.0, 1.0, 1.0};
    boundary_sources[0] = make_shared<Boundary_Source>(index,
                                                       dependencies,
                                                       angular,
                                                       energy,
                                                       boundary_source,
                                                       alpha);

    return boundary_sources;
}

void
get_solid_geometry(vector<shared_ptr<Material> > materials,
                   vector<shared_ptr<Boundary_Source> > boundary_sources,
                   shared_ptr<Constructive_Solid_Geometry> &solid,
                   vector<shared_ptr<Cartesian_Plane> > &boundary_surfaces)
{
    // Initialize surfaces
    int dimension = 2;
    int number_of_surfaces = 5;
    int number_of_boundary_surfaces = 4;
    vector<shared_ptr<Surface> > surfaces(number_of_surfaces);
    boundary_surfaces.resize(number_of_boundary_surfaces);
    
    // Get boundary planes
    {
        vector<int> surface_dimensions = {0, 0, 1, 1};
        vector<double> surface_positions = {-1.2, 1.2, -1.2, 1.2};
        vector<double> surface_normals = {-1.0, 1.0, -1.0, 1.0};
        for (int i = 0; i < 4; ++i)
        {
            boundary_surfaces[i]
                = make_shared<Cartesian_Plane>(i, // index
                                               dimension,
                                               Surface::Surface_Type::BOUNDARY,
                                               surface_dimensions[i],
                                               surface_positions[i],
                                               surface_normals[i]);
            surfaces[i] = boundary_surfaces[i];
            surfaces[i]->set_boundary_source(boundary_sources[0]);
        }
    }

    // Get internal cylinder
    {
        int i = 4;
        double radius = 0.5;
        vector<double> origin = {0.0, 0.0};
        surfaces[i] = make_shared<Cylinder_2D>(i,
                                               Surface::Surface_Type::INTERNAL,
                                               radius, 
                                               origin); 
    }

    // Initialize regions
    int number_of_regions = 2;
    vector<shared_ptr<Region> > regions(number_of_regions);

    // Get pincell region
    {
        int i = 0;
        vector<Surface::Relation> surface_relations
            = {Surface::Relation::INSIDE};
        vector<shared_ptr<Surface> > local_surfaces
            = {surfaces[4]};
        
        regions[i] = make_shared<Region>(i,
                                         materials[i],
                                         surface_relations,
                                         local_surfaces);
    }
    
    // Get moderator region
    {
        int i = 1;
        vector<Surface::Relation> surface_relations
            = {Surface::Relation::NEGATIVE,
               Surface::Relation::NEGATIVE,
               Surface::Relation::NEGATIVE,
               Surface::Relation::NEGATIVE,
               Surface::Relation::OUTSIDE};
        
        regions[i] = make_shared<Region>(i, // index
                                         materials[i],
                                         surface_relations,
                                         surfaces);
    }

    // Get solid geometry
    solid =  make_shared<Constructive_Solid_Geometry>(dimension,
                                                      surfaces,
                                                      regions,
                                                      materials,
                                                      boundary_sources);
}

shared_ptr<Weak_Spatial_Discretization_Options>
get_weak_options(bool supg,
                 bool use_flux,
                 int num_dimensional_points,
                 shared_ptr<Angular_Discretization> angular,
                 shared_ptr<Energy_Discretization> energy,
                 shared_ptr<Solid_Geometry> solid)
{
    shared_ptr<Weak_Spatial_Discretization_Options> weak_options
        = make_shared<Weak_Spatial_Discretization_Options>();

    weak_options->integration_ordinates = 16;
    weak_options->limits = {{-1.2, 1.2}, {-1.2, 1.2}};
    weak_options->solid = solid;
    weak_options->dimensional_cells = {2 * (num_dimensional_points - 1), 2 * (num_dimensional_points - 1)};
    weak_options->identical_basis_functions = Weak_Spatial_Discretization_Options::Identical_Basis_Functions::TRUE;
    
    if (supg)
    {
        weak_options->include_supg = true;
    }

    if (use_flux)
    {
        int number_of_groups = energy->number_of_groups();
        int number_of_moments = angular->number_of_moments();
        int number_of_points = num_dimensional_points * num_dimensional_points;
        weak_options->weighting = Weak_Spatial_Discretization_Options::Weighting::FLUX;
        weak_options->flux_coefficients.assign(number_of_groups * number_of_moments * number_of_points, 1.);
    }

    return weak_options;
}

void
get_discretizations(bool supg,
                    bool use_flux,
                    int num_dimensional_points,
                    double tau,
                    shared_ptr<Angular_Discretization> &angular,
                    shared_ptr<Energy_Discretization> &energy,
                    shared_ptr<Constructive_Solid_Geometry> &solid,
                    shared_ptr<Weak_Spatial_Discretization> &spatial)
{
    // Initialize angular discretization
    int dimension = 2;
    int number_of_scattering_moments = 2;
    angular = make_shared<LDFE_Quadrature>(dimension,
                                           number_of_scattering_moments,
                                           2); // rule

    // Initialize energy discretization
    int number_of_groups = 3;
    energy = make_shared<Energy_Discretization>(number_of_groups);

    // Initialize materials
    vector<shared_ptr<Material> > materials = get_materials(angular,
                                                            energy);
    
    // Initialize boundary sources
    vector<shared_ptr<Boundary_Source> > boundary_sources = get_boundary_sources(angular,
                                                                                 energy);

    // Initialize solid geometry
    vector<shared_ptr<Cartesian_Plane> > boundary_surfaces;
    get_solid_geometry(materials,
                       boundary_sources,
                       solid,
                       boundary_surfaces);
    
    // Initialize weight function options
    shared_ptr<Weight_Function_Options> weight_options
        = make_shared<Weight_Function_Options>();
    weight_options->tau_const = tau;
    
    // Initialize spatial options
    shared_ptr<Weak_Spatial_Discretization_Options> weak_options
        = get_weak_options(supg,
                           use_flux,
                           num_dimensional_points,
                           angular,
                           energy,
                           solid);

    // Initialize spatial discretization
    Weak_Spatial_Discretization_Factory spatial_factory(solid,
                                                        boundary_surfaces);
    
    spatial
        = spatial_factory.get_simple_discretization(num_dimensional_points, // num dimensional points
                                                    3, // radius num intervals
                                                    true, // basis mls
                                                    true, // weight mls
                                                    "wendland11", // basis type
                                                    "wendland11", // weight type
                                                    weight_options,
                                                    weak_options);
}

shared_ptr<Vector_Operator>
get_standard_operator(bool use_flux,
                      int num_dimensional_points)
{
    // Get discretizations
    shared_ptr<Angular_Discretization> angular;
    shared_ptr<Energy_Discretization> energy;
    shared_ptr<Constructive_Solid_Geometry> solid;
    shared_ptr<Weak_Spatial_Discretization> spatial;
    
    get_discretizations(false, // supg
                        use_flux,
                        num_dimensional_points,
                        0.0, // tau
                        angular,
                        energy,
                        solid,
                        spatial);

    // Get individual operators
    shared_ptr<Moment_To_Discrete> M
        = make_shared<Moment_To_Discrete>(spatial,
                                          angular,
                                          energy);
    shared_ptr<Scattering> S
        = make_shared<Scattering>(spatial,
                                  angular,
                                  energy);
    shared_ptr<Fission> F
        = make_shared<Fission>(spatial,
                               angular,
                               energy);
    shared_ptr<Moment_Weighting_Operator> W
        = make_shared<Moment_Weighting_Operator>(spatial,
                                                 angular,
                                                 energy);

    // Combine combined operator
    return M * (S + F) * W;
}
                                                  
shared_ptr<Vector_Operator>
get_supg_operator(bool use_flux,
                  int num_dimensional_points,
                  double tau)
{
    // Get discretizations
    shared_ptr<Angular_Discretization> angular;
    shared_ptr<Energy_Discretization> energy;
    shared_ptr<Constructive_Solid_Geometry> solid;
    shared_ptr<Weak_Spatial_Discretization> spatial;
    
    get_discretizations(true, // supg
                        use_flux,
                        num_dimensional_points,
                        tau,
                        angular,
                        energy,
                        solid,
                        spatial);

    // Return different operators depending on flux option
    if (use_flux)
    {
        // Get individual operators
        shared_ptr<Vector_Operator> S
            = make_shared<Combined_SUPG_Scattering>(spatial,
                                                    angular,
                                                    energy);
        shared_ptr<Vector_Operator> F
            = make_shared<Combined_SUPG_Fission>(spatial,
                                                 angular,
                                                 energy);
        shared_ptr<Vector_Operator> W
            = make_shared<Moment_Weighting_Operator>(spatial,
                                                     angular,
                                                     energy);
        // Return combined operator
        return (S + F) * W;
    }
    else
    {
        // Get individual operators
        shared_ptr<Vector_Operator> M2
            = make_shared<SUPG_Moment_To_Discrete>(spatial,
                                                   angular,
                                                   energy,
                                                   true); // include double dimensional moments
        shared_ptr<Vector_Operator> S
            = make_shared<SUPG_Scattering>(spatial,
                                           angular,
                                           energy);
        shared_ptr<Vector_Operator> F
            = make_shared<SUPG_Fission>(spatial,
                                        angular,
                                        energy);
        shared_ptr<Vector_Operator> W
            = make_shared<Moment_Weighting_Operator>(spatial,
                                                     angular,
                                                     energy);
        shared_ptr<Vector_Operator> N
            = make_shared<Discrete_Normalization_Operator>(spatial,
                                                           angular,
                                                           energy);
        
        // Return combined operator
        return N * M2 * (S + F) * W;
    }
}

int main()
{
    int checksum = 0;

    
    
    return checksum;
}
