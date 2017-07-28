#include <iomanip>
#include <iostream>
#include <vector>

#include "Angular_Discretization_Factory.hh"
#include "Boundary_Source.hh"
#include "Cartesian_Plane.hh"
#include "Check_Equality.hh"
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
namespace ce = Check_Equality;

// Create a three-group, anisotropically-scattering pincell problem with reflective boundaries
// Test whether the various scattering and fission operators return identical results when they should

Random_Number_Generator<double> rng(0, // lower bound
                                    1, // upper bound
                                    203); // seed

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
                                           1); // rule

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

void
get_standard_operator(bool use_flux,
                      int num_dimensional_points,
                      shared_ptr<Angular_Discretization> &angular,
                      shared_ptr<Energy_Discretization> &energy,
                      shared_ptr<Constructive_Solid_Geometry> &solid,
                      shared_ptr<Weak_Spatial_Discretization> &spatial,
                      shared_ptr<Vector_Operator> &oper)
{
    // Get discretizations
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
    oper = M * (S + F) * W;
}
                                                  
void
get_supg_operator(bool use_flux,
                  int num_dimensional_points,
                  double tau,
                  shared_ptr<Angular_Discretization> &angular,
                  shared_ptr<Energy_Discretization> &energy,
                  shared_ptr<Constructive_Solid_Geometry> &solid,
                  shared_ptr<Weak_Spatial_Discretization> &spatial,
                  shared_ptr<Vector_Operator> &oper)
{
    // Get discretizations
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
        oper = (S + F) * W;
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
        oper = N * M2 * (S + F) * W;
    }
}

void
get_random_coefficients(shared_ptr<Angular_Discretization> angular,
                        shared_ptr<Energy_Discretization> energy,
                        shared_ptr<Spatial_Discretization> spatial,
                        vector<double> &coefficients)
{
    // Set size of coefficient vector
    int number_of_moments = angular->number_of_moments();
    int number_of_groups = energy->number_of_groups();
    int number_of_points = spatial->number_of_points();
    coefficients.resize(number_of_moments * number_of_groups * number_of_points);
    
    for (int i = 0; i < number_of_points; ++i)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            for (int m = 0; m < number_of_moments; ++m)
            {
                // Get random number
                double rand = rng.scalar();

                // Scale random number for non-scalar coefficients
                if (m != 0)
                {
                    // Allow random number to be negative
                    rand = 2 * rand - 1;

                    // Ensure that random number is not larger than a third of the scalar coefficient
                    int k0 = g + number_of_groups * (0 + number_of_moments * i);
                    rand *= 0.3 * coefficients[k0];
                }
                
                int k = g + number_of_groups * (m + number_of_moments * i);
                coefficients[k] = rand;
            }
        }
    }
}

void
print_error(shared_ptr<Angular_Discretization> angular,
            shared_ptr<Energy_Discretization> energy,
            shared_ptr<Spatial_Discretization> spatial,
            double tolerance,
            string desc1,
            vector<double> val1,
            string desc2,
            vector<double> val2)
{
    // Check size data
    int number_of_points = spatial->number_of_points();
    int number_of_groups = energy->number_of_groups();
    int number_of_ordinates = angular->number_of_ordinates();
    int flux_size = number_of_points * number_of_groups * number_of_ordinates;
    Assert(val1.size() == flux_size);
    Assert(val2.size() == flux_size);
    
    // Print header
    int w1 = 10;
    int w2 = 16;
    cout << setw(w1) << "cell";
    cout << setw(w1) << "group";
    cout << setw(w1) << "ordinate";
    cout << setw(w2) << desc1;
    cout << setw(w2) << desc2;
    cout << setw(w2) << "difference";
    cout << endl;

    // Check points for error
    for (int i = 0; i < number_of_points; ++ i)
    {
        for (int g = 0; g < number_of_groups; ++g)
        {
            for (int o = 0; o < number_of_ordinates; ++o)
            {
                // Get difference between values
                int k = g + number_of_groups * (o + number_of_ordinates * i);
                double diff = val1[k] - val2[k];

                // If difference is too large, print data
                if (abs(diff) > tolerance)
                {
                    cout << setw(w1) << i;
                    cout << setw(w1) << g;
                    cout << setw(w1) << o;
                    cout << setw(w2) << val1[k];
                    cout << setw(w2) << val2[k];
                    cout << setw(w2) << diff;
                    cout << endl;
                }
            }
        }
    }
}
    
int check_supg_operators(int num_dimensional_points)
{
    // Set preliminary values
    int checksum = 0;
    double tau = 0.9;
    cout << "check_supg_operators running for ";
    cout << num_dimensional_points;
    cout << " dimensional points";
    cout << endl;
        
    
    // Get simple SUPG operator
    cout << "creating simple SUPG operator" << endl;
    shared_ptr<Angular_Discretization> simple_angular;
    shared_ptr<Energy_Discretization> simple_energy;
    shared_ptr<Constructive_Solid_Geometry> simple_solid;
    shared_ptr<Weak_Spatial_Discretization> simple_spatial;
    shared_ptr<Vector_Operator> simple_oper;
    get_supg_operator(false, // use_flux
                      num_dimensional_points,
                      tau,
                      simple_angular,
                      simple_energy,
                      simple_solid,
                      simple_spatial,
                      simple_oper);

    // Get flux-weighted SUPG operator
    cout << "creating flux-weighted SUPG operator" << endl;
    shared_ptr<Angular_Discretization> flux_angular;
    shared_ptr<Energy_Discretization> flux_energy;
    shared_ptr<Constructive_Solid_Geometry> flux_solid;
    shared_ptr<Weak_Spatial_Discretization> flux_spatial;
    shared_ptr<Vector_Operator> flux_oper;
    get_supg_operator(true, // use_flux
                      num_dimensional_points,
                      tau,
                      flux_angular,
                      flux_energy,
                      flux_solid,
                      flux_spatial,
                      flux_oper);

    // Ensure that size information is the same for the two discretizations
    Assert(simple_angular->number_of_moments() == flux_angular->number_of_moments());
    Assert(simple_energy->number_of_groups() == flux_energy->number_of_groups());
    Assert(simple_spatial->number_of_points() == flux_spatial->number_of_points());
    
    // Get random coefficients
    vector<double> coefficients;
    get_random_coefficients(simple_angular,
                            simple_energy,
                            simple_spatial,
                            coefficients);

    // Apply the operators to the random coefficients
    cout << "applying operators" << endl;
    vector<double> simple_result = coefficients;
    (*simple_oper)(simple_result);
    vector<double> flux_result = coefficients;
    (*flux_oper)(flux_result);

    // Check to ensure that the two operators returned the same results
    double tolerance = 1e-10;
    if (ce::approx(simple_result, flux_result, tolerance))
    {
        cout << "test_passed" << endl;
    }
    else
    {
        cout << "supg operator results differ" << endl;
        print_error(simple_angular,
                    simple_energy,
                    simple_spatial,
                    tolerance,
                    "simple",
                    simple_result,
                    "flux",
                    flux_result);
        checksum += 1;
    }

    return checksum;
}

int main()
{
    int checksum = 0;
    
    checksum += check_supg_operators(11);
    
    return checksum;
}
