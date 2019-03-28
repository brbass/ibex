#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

#include "Heat_Transfer_Integration.hh"

using namespace std;

int test_integral_2d()
{
    // Set data
    double length1 = 1.;
    double length2 = 2.;
    double conduction1 = 0.01;
    double conduction2 = 0.01;
    double convection = 10.;
    double source1 = 10.;
    double source2 = 10.;
    double temperature_inf = 600.;
    int dimension = 2;
    int num_points = 2;
    
    // Get solid geometry
    int dimension = 2;
    vector<vector<double> > limits = {{-length2, length2},
                                      {-length2, length2}};
    shared_ptr<Solid_Geometry> solid;
    vector<shared_ptr<Cartesian_Plane> > surfaces;
    Heat_Transfer_Factory factory;
    factory.get_solid(dimension,
                      limits,
                      solid,
                      surfaces);
    
    // Get data
    shared_ptr<Constant_Heat_Transfer_Data> data
        = make_shared<Constant_Heat_Transfer_Data>(dimension,
                                                   length1,
                                                   length2,
                                                   conduction1,
                                                   conduction2,
                                                   convection,
                                                   source1,
                                                   source2,
                                                   temperature_inf);

    // Get weak spatial discretization options
    shared_ptr<Dimensional_Moments> dimensional_moments
        = make_shared<Dimensional_Moments>(false, // supg
                                           dimension);
    shared_ptr<Weight_Function_Options> weight_options
        = make_shared<Weight_Function_Options>();
    shared_ptr<Weak_Spatial_Discretization_Options> options
        = make_shared<Weak_Spatial_Discretization_Options>();
    options->weighting = Weak_Spatial_Discretization_Options::Weighting::FLAT;
    options->external_integral_calculation = false;
    options->include_supg = false;
    options->tau_scaling = Weak_Spatial_Discretization_Options::Tau_Scaling::NONE;
    options->identical_basis_functions = Weak_Spatial_Discretization_Options::Identical_Basis_Functions::TRUE;
    options->perform_integration = false;
    options->adaptive_quadrature = false;

    // Get meshless functions
    vector<double> radii = {3.0, 4.0};
    vector<vector<double> > points = {{0.0, 0.0}, {0.5, 0.3}};
    RBF_Factory rbf_factory;
    shared_ptr<RBF> rbf = rbf_factory.get_rbf("wendland11");
    shared_ptr<Distance> distance = make_shared<Cartesian_Distance>(dimension);
    Meshless_Function_Factory meshless_factory;
    vector<shared_ptr<Meshless_Function> > functions;
    meshless_factory.get_rbf_functions(num_points,
                                       radii,
                                       points,
                                       rbf,
                                       distance,
                                       functions);
    
    // Get bases and weights
    vector<shared_ptr<Basis_Function> > bases(num_points);
    vector<shared_ptr<Weight_Function> > weights(num_points);
    for (int i = 0; i < num_points; ++i)
    {
        shared_ptr<Meshless_Function> meshless_function
            = make_shared<
        bases[i] = make_shared<Basis_Function>(i,
                                               dimension,
                                               meshless_function,
                                               surfaces);
        weights[i] = make_shared<Weight_Function>(i,
                                                  dimension,
                                                  weight_options,
                                                  weak_options,
                                                  meshless_function,
                                                  dimensional_moments,
                                                  solid_geometry,
                                                  surfaces);
    }
    
    // Get spatial discretization
    shared_ptr<Weak_Spatial_Discretization> spatial
        = make_shared<Weak_Spatial_Discretization>(bases,
                                                   weights,
                                                   dimensional_moments,
                                                   options);
    
    // Get heat transfer integration options
    shared_ptr<Heat_Transfer_Integration_Options> integration_options
        = make_shared<Heat_Transfer_Integration_Options>();
    integration_options->geometry = Heat_Transfer_Integration_Options::Geometry::CYLINDRICAL_2D;
    
    // Perform integration
    shared_ptr<Heat_Transfer_Integration> integration
        = make_shared<Heat_Transfer_Integration>(integration_options,
                                                 data,
                                                 spatial);
    vector<vector<double> > matrix = integration->matrix();
    for (int i = 0; i < num_points; ++i)
    {
        for (int j = 0; j < num_points; ++j)
        {
            cout << matrix[i][j] << "\t";
        }
        cout << endl;
    }
}


int main(int argc, char **argv)
{
    int checksum = 0;
    
    MPI_Init(&argc, &argv);
    
    test_integral_2d();
    
    MPI_Finalize();
    
    return checksum;
}
