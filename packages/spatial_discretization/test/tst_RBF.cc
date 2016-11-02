#include <iostream>
#include <memory>
#include <vector>

#include "Check_Equality.hh"
#include "Gaussian_RBF.hh"
#include "Inverse_Multiquadric_RBF.hh"
#include "Multiquadric_RBF.hh"
#include "RBF.hh"
#include "Wendland_RBF.hh"

namespace ce = Check_Equality;

using namespace std;

int test_rbf(shared_ptr<RBF> rbf,
             int number_of_cases,
             string description,
             vector<double> const &points,
             vector<double> const &expected_basis,
             vector<double> const &expected_d_basis,
             vector<double> const &expected_dd_basis)
{
    int checksum = 0;
    
    for (int i = 0; i < number_of_cases; ++i)
    {
        if (!ce::approx(rbf->basis(points[i]), expected_basis[i], 1e-15))
        {
            cout << description << " basis failed at " << points[i] << endl;
            checksum += 1;
        }

        if (!ce::approx(rbf->d_basis(points[i]), expected_d_basis[i], 1e-15))
        {
            cout << description << " d_basis failed at " << points[i] << endl;
            checksum += 1;
        }

        if (!ce::approx(rbf->dd_basis(points[i]), expected_dd_basis[i], 1e-15))
        {
            cout << description << " dd_basis failed at " << points[i] << endl;
            checksum += 1;
        }
    }

    return checksum;
}

int main()
{
    int checksum = 0;

    int const number_of_cases = 3;
    vector<double> const points = {0., 0.5, 1.0};
    
    checksum += test_rbf(make_shared<Gaussian_RBF>(),
                         number_of_cases,
                         "gaussian",
                         points,
                         {1., exp(-1. / 4.), exp(-1.)},
                         {0., -exp(-1. / 4.), -2 * exp(-1.)},
                         {-2, -exp(-1. / 4.), 2 * exp(-1.)});
    checksum += test_rbf(make_shared<Inverse_Multiquadric_RBF>(),
                         number_of_cases,
                         "inverse_multiquadric",
                         points,
                         {1., 2. / sqrt(5.), 1. / sqrt(2.)},
                         {0, -4. / (5. * sqrt(5.)), -1. / (2. * sqrt(2.))},
                         {-1., -16. / (25. * sqrt(5.)), 1. / (4. * sqrt(2.))});
    checksum += test_rbf(make_shared<Multiquadric_RBF>(),
                         number_of_cases,
                         "multiquadric",
                         points,
                         {1., sqrt(5.) / 2., sqrt(2.)},
                         {0, 1. / sqrt(5.), 1. / sqrt(2.)},
                         {1, 8. / (5. * sqrt(5.)), 1. / (2. * sqrt(2.))});
    checksum += test_rbf(make_shared<Wendland_RBF>(0),
                         number_of_cases,
                         "wendland0",
                         points,
                         {1., 1. / 4., 0.},
                         {-2., -1., 0.},
                         {2., 2., 2.});
    checksum += test_rbf(make_shared<Wendland_RBF>(1),
                         number_of_cases,
                         "wendland1",
                         points,
                         {1., 3. / 16., 0.},
                         {0., -5. / 4., 0.},
                         {-20., 5., 0.});
    checksum += test_rbf(make_shared<Wendland_RBF>(2),
                         number_of_cases,
                         "wendland2",
                         points,
                         {3., 83. / 256., 0.},
                         {0., -49. / 16., 0.},
                         {-56., 161. / 8., 0.});
    checksum += test_rbf(make_shared<Wendland_RBF>(3),
                         number_of_cases,
                         "wendland3",
                         points,
                         {1., 61. / 1024., 0.},
                         {0., -187. / 256., 0.},
                         {-22., 869. / 128., 0.});
        
    return checksum;
}
