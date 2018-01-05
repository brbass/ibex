#include <iostream> 
#include <limits>

#include "Check_Equality.hh"
#include "Cross_Section.hh"
#include "Discrete_To_Moment.hh"
#include "Energy_Discretization.hh"
#include "Gauss_Legendre_Quadrature.hh"
#include "LDFE_Quadrature.hh"
#include "Material.hh"
#include "Material_Factory.hh"
#include "Moment_To_Discrete.hh"
#include "Simple_Point.hh"
#include "Random_Number_Generator.hh"
#include "Simple_Spatial_Discretization.hh"
#include "Vector_Operator_Functions.hh"

namespace ce = Check_Equality;

using namespace std;

// Get spatial, angular and energy discretizations
void get_discretization(int dimension,
                        int number_of_points,
                        int number_of_groups,
                        int quadrature_rule,
                        int number_of_scattering_moments,
                        shared_ptr<Spatial_Discretization> &spatial,
                        shared_ptr<Angular_Discretization> &angular,
                        shared_ptr<Energy_Discretization> &energy)
{
    // Get energy discretization
    energy = make_shared<Energy_Discretization>(number_of_groups);

    // Get angular discretization
    switch(dimension)
    {
    case 1:
        angular = make_shared<Gauss_Legendre_Quadrature>(dimension,
                                                         number_of_scattering_moments,
            quadrature_rule);
        break;
    default:
        angular = make_shared<LDFE_Quadrature>(dimension,
                                               number_of_scattering_moments,
                                               quadrature_rule);
        break;
    }
    
    // Create dummy material
    Cross_Section::Dependencies none_group;
    none_group.energy = Cross_Section::Dependencies::Energy::GROUP;
    Cross_Section::Dependencies scattering_group2;
    scattering_group2.angular = Cross_Section::Dependencies::Angular::SCATTERING_MOMENTS;
    scattering_group2.energy = Cross_Section::Dependencies::Energy::GROUP_TO_GROUP;

    vector<double> sigma_t_data(number_of_groups, 1);
    vector<double> sigma_s_data(number_of_groups * number_of_groups * number_of_scattering_moments, 0);
    vector<double> nu_data(number_of_groups, 0);
    vector<double> sigma_f_data(number_of_groups, 0);
    vector<double> chi_data(number_of_groups, 0);
    vector<double> internal_source_data(number_of_groups, 0);

    Material_Factory factory(angular,
                             energy);
    
    shared_ptr<Material> material
        = factory.get_standard_material(0,
                                        sigma_t_data,
                                        sigma_s_data,
                                        nu_data,
                                        sigma_f_data,
                                        chi_data,
                                        internal_source_data);
    
    // Create dummy spatial discretization    
    vector<shared_ptr<Point> > points(number_of_points);

    for (int i = 0; i < number_of_points; ++i)
    {
        vector<double> position(dimension, 0);
        points[i] = make_shared<Simple_Point>(i,
                                              dimension,
                                              Point::Point_Type::INTERNAL,
                                              material,
                                              position);
    }
    
    spatial = make_shared<Simple_Spatial_Discretization>(points);
}

// Can convert from phi to psi to phi without losing information
bool lossless_phi(int dimension,
                  int quadrature_rule,
                  int number_of_scattering_moments)
{
    switch(dimension)
    {
    case 1:
        return 2 * number_of_scattering_moments - 1 <= quadrature_rule;
    default:
        switch(quadrature_rule)
        {
        case 1:
            return number_of_scattering_moments <= quadrature_rule;
        case 2:
            return number_of_scattering_moments <= quadrature_rule;
        default:
            return number_of_scattering_moments < 3;
        }
    }
}

// Test a single case of moment-to-discrete and reverse
int test_moment_discrete(int dimension,
                         int quadrature_rule,
                         int number_of_scattering_moments)
{
    int checksum = 0;
    
    double const tolerance = 1000 * numeric_limits<double>::epsilon();

    // Initialize data
    int number_of_groups = 2;
    int number_of_points = 10;
    
    shared_ptr<Spatial_Discretization> spatial;
    shared_ptr<Angular_Discretization> angular;
    shared_ptr<Energy_Discretization> energy;
    
    get_discretization(dimension,
                       number_of_points,
                       number_of_groups,
                       quadrature_rule,
                       number_of_scattering_moments,
                       spatial,
                       angular,
                       energy);

    // Get moment to discrete and discrete to moment operators
    int number_of_moments = angular->number_of_moments();
    int number_of_ordinates = angular->number_of_ordinates();
    int phi_size = number_of_points * number_of_groups * number_of_moments;
    
    shared_ptr<Vector_Operator> M
        = make_shared<Moment_To_Discrete>(spatial,
                                          angular,
                                          energy);
    shared_ptr<Vector_Operator> D
        = make_shared<Discrete_To_Moment>(spatial,
                                          angular,
                                          energy);

    // Get combined operators for testing
    shared_ptr<Vector_Operator> MD = M * D;
    shared_ptr<Vector_Operator> DM = D * M;
    
    Random_Number_Generator<double> rng(-1, 1, 9834);

    int number_of_tests = 10;

    // Test for a number of random cases
    for (int i = 0; i < number_of_tests; ++i)
    {
        // Get random moment data
        vector<double> const phi = rng.vector(phi_size);
        vector<double> phi_new(phi);

        // Convert from phi to psi and back to phi
        (*DM)(phi_new);

        // Check whether original and final are the same
        if (!ce::approx(phi, phi_new, tolerance))
        {
            cout << "moment discrete from phi failed for test ";
            cout << i;
            cout << endl;
            cout << "\tdimension: ";
            cout << dimension;
            cout << "\tmoments: ";
            cout << number_of_scattering_moments;
            cout << "\tquadrature_rule: ";
            cout << quadrature_rule;
            cout << endl;
            checksum += 1;
        }
    }
    
    return checksum;
}

// Run all tests
int main()
{
    int checksum = 0;

    // Check one dimension
    {
        int dimension = 1;
        
        for (int moments = 1; moments < 10; ++moments)
        {
            for (int quad_rule = 2; quad_rule < 20; ++++quad_rule)
            {
                // Check whether the conversion can be done without loss
                if (lossless_phi(dimension,
                                 quad_rule,
                                 moments))
                {
                    checksum += test_moment_discrete(dimension,
                                                     quad_rule,
                                                     moments);
                }
            }
        }
        
    }

    // Check two and three dimensions
    for (int dimension = 2; dimension <= 3; ++dimension)
    {
        for (int quad_rule = 1; quad_rule < 7; ++quad_rule)
        {
            for (int moments = 1; moments < 4; ++moments)
            {
                // Check whether the conversion can be done without loss
                if (lossless_phi(dimension,
                                 quad_rule,
                                 moments))
                {
                    checksum += test_moment_discrete(dimension,
                                                     quad_rule,
                                                     moments);
                }
            }
        }
    }
        
    return checksum;
}
