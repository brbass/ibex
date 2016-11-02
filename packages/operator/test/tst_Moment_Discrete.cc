#include <iostream> 
#include <limits>

#include "Check_Equality.hh"
#include "Discrete_To_Moment.hh"
#include "Energy_Discretization.hh"
#include "Gauss_Legendre_Quadrature.hh"
#include "LDFE_Quadrature.hh"
#include "Moment_To_Discrete.hh"
#include "Point.hh"
#include "Random_Number_Generator.hh"
#include "Simple_Spatial_Discretization.hh"
#include "Vector_Operator_Functions.hh"

namespace ce = Check_Equality;

using namespace std;

void get_discretization(int dimension,
                        int number_of_points,
                        int number_of_groups,
                        int quadrature_rule,
                        int number_of_scattering_moments,
                        shared_ptr<Spatial_Discretization> &spatial,
                        shared_ptr<Angular_Discretization> &angular,
                        shared_ptr<Energy_Discretization> &energy)
{
    // Energy discretization
    
    energy = make_shared<Energy_Discretization>(number_of_groups);

    // Angular discretization
    
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

    // Material

    vector<double> sigma_t(number_of_groups, 1);
    vector<double> sigma_s(number_of_groups * number_of_groups * number_of_scattering_moments, 0);
    vector<double> nu(number_of_groups, 0);
    vector<double> sigma_f(number_of_groups, 0);
    vector<double> chi(number_of_groups, 0);
    vector<double> internal_source(number_of_groups, 0);
    
    shared_ptr<Material> material
        = make_shared<Material>(0,
                                angular,
                                energy,
                                sigma_t,
                                sigma_s,
                                nu,
                                sigma_f,
                                chi,
                                internal_source);
    
    // Spatial discretization
    
    vector<shared_ptr<Point> > points(number_of_points);

    for (int i = 0; i < number_of_points; ++i)
    {
        vector<double> position(dimension, 0);
        points[i] = make_shared<Point>(i,
                                       dimension,
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
            return number_of_scattering_moments < quadrature_rule;
        }
    }
}


int test_moment_discrete(int dimension,
                         int quadrature_rule,
                         int number_of_scattering_moments)
{
    int checksum = 0;

    double const tolerance = 1000 * numeric_limits<double>::epsilon();
    
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
    
    shared_ptr<Vector_Operator> MD = M * D;
    shared_ptr<Vector_Operator> DM = D * M;

    Random_Number_Generator<double> rng(-1, 1, 9834);

    int number_of_tests = 10;
    
    for (int i = 0; i < number_of_tests; ++i)
    {
        vector<double> const phi = rng.vector(phi_size);
        vector<double> phi_new(phi);
            
        (*DM)(phi_new);
            
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

int main()
{
    int checksum = 0;

    {
        int dimension = 1;
        
        for (int moments = 1; moments < 10; ++moments)
        {
            for (int quad_rule = 2; quad_rule < 20; ++++quad_rule)
            {
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

    for (int dimension = 2; dimension < 4; ++dimension)
    {
        for (int quad_rule = 1; quad_rule < 4; ++quad_rule)
        {
            for (int moments = 1; moments < 4; ++moments)
            {
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
