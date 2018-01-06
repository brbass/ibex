#include <cmath>
#include <iomanip>
#include <iostream>

#include "Angular_Discretization.hh"
#include "Angular_Discretization_Factory.hh"
#include "Check_Equality.hh"

namespace ce = Check_Equality;

using namespace std;

int check_coefficients(int dimension,
                       int number_of_scattering_moments,
                       int angular_rule,
                       vector<vector<vector<double> > > const &expected_coefficients)
{
    int checksum = 0;
    
    // Get angular discretization
    Angular_Discretization_Factory factory;
    shared_ptr<Angular_Discretization> angular
        = factory.get_angular_discretization(dimension,
                                             number_of_scattering_moments,
                                             angular_rule);
    int number_of_moments = angular->number_of_moments();
    
    // Get coefficients
    vector<vector<int> > indices;
    vector<vector<vector<double> > > coefficients;
    angular->manufactured_coefficients(indices,
                                       coefficients);

    // Print
    vector<int> l_vals = angular->harmonic_degrees();
    vector<int> m_vals = angular->harmonic_orders();
    for (int o1 = 0; o1 < number_of_moments; ++o1)
    {
        int l1 = l_vals[o1];
        int m1 = m_vals[o1];

        int num_indices = indices[o1].size();
        cout << o1 << endl;
        for (int n = 0; n < num_indices; ++n)
        {
            int o2 = indices[o1][n];
            int l2 = l_vals[o2];
            int m2 = m_vals[o2];

            int w = 4;
            cout << setw(w) << l1;
            cout << setw(w) << m1;
            cout << setw(w) << l2;
            cout << setw(w) << m2;
            for (int d = 0; d < dimension; ++d)
            {
                cout << setw(12) << coefficients[o1][n][d];
            }
            cout << endl;
        }
    }
    cout << endl;

    return checksum;
}

int main()
{
    int checksum = 0;

    int dimension;
    int angular_rule;
    int scattering_order;
    vector<vector<vector<double> > > coefficients;
    {
        dimension = 1;
        angular_rule = 16;
        scattering_order = 4;
        // coefficients = {{{}}};
        checksum += check_coefficients(dimension,
                                       scattering_order,
                                       angular_rule,
                                       coefficients);
    }
    {
        dimension = 3;
        angular_rule = 4;
        scattering_order = 3;

        double n13 = 1./3.;
        double n23 = 2./3.;
        double s13 = 1./sqrt(3.);
        // coefficients = {{{0, 0, 1},
        //                  {1, 0, 0},
        //                  {0, 1, 0}},
        //                 {{0, 0, n13},
        //                  {0, s13, 0},
        //                  {s13, 0, 0},
        //                  {0, 0, -third},
        //                  {0., third, -s13}},
        //                 {{third, 0, 0},
        //                  {0, 0, s13},
        //                  {n23, 0, 0},
        //                  {0, s13, 0}},
        //                 {{0, n13, 0},
        //                  {};
        checksum += check_coefficients(dimension,
                                       scattering_order,
                                       angular_rule,
                                       coefficients);
    }

    return checksum;
}
