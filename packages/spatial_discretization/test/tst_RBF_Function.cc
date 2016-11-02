#include <cmath>
#include <iostream>
#include <memory>

#include "Cartesian_Distance.hh"
#include "Check_Equality.hh"
#include "Derivative_RBF_Function.hh"
#include "Distance.hh"
#include "Multiquadric_RBF.hh"
#include "RBF.hh"
#include "RBF_Function.hh"
#include "XML_Functions.hh"

using namespace std;

namespace ce = Check_Equality;
namespace xf = XML_Functions;

int test_rbf_function(shared_ptr<RBF_Function> rbf_function,
                      string test_case,
                      int dimension,
                      double const shape,
                      double const expected_basis,
                      vector<double> const &expected_grad,
                      vector<double> const &expected_double_grad,
                      vector<double> const &r,
                      vector<double> const &r0,
                      vector<double> const &direction)
{
    int checksum = 0;

    int group = 0;
    double tol = 1e-15;
    
    // Check basis
    
    double basis = rbf_function->basis(group,
                                       shape,
                                       r,
                                       r0,
                                       direction);

    if (!ce::approx(basis, expected_basis, tol))
    {
        cout << "basis failed for ";
        cout << test_case;
        cout << endl;
        cout << "\texpected: ";
        cout << expected_basis;
        cout << "\tcalculated: ";
        cout << basis;
        cout << endl;
        checksum += 1;
    }

    // Check first derivatives
    
    if (rbf_function->derivative_available(1))
    {
        for (int d = 0; d < dimension; ++d)
        {
            double d_basis = rbf_function->d_basis(group,
                                                   d,
                                                   shape,
                                                   r,
                                                   r0,
                                                   direction);
            
            if (!ce::approx(d_basis, expected_grad[d], tol))
            {
                cout << "d_basis in dimension ";
                cout << d;
                cout << " failed for ";
                cout << test_case;
                cout << endl;
                cout << "\texpected: ";
                cout << expected_grad[d];
                cout << "\tcalculated: ";
                cout << d_basis;
                cout << endl;
                checksum += 1;
            }
        }
        
        vector<double> grad = rbf_function->gradient_basis(group,
                                                           shape,
                                                           r,
                                                           r0,
                                                           direction);

        if (!ce::approx(grad, expected_grad, tol))
        {
            string eg, gr;
            xf::vector_to_string(eg, expected_grad);
            xf::vector_to_string(gr, grad);
            
            cout << "grad basis failed for ";
            cout << test_case;
            cout << endl;
            cout << "\texpected: ";
            cout << eg;
            cout << endl;
            cout << "\tcalculated: ";
            cout << gr;
            cout << endl;
            checksum += 1;
        }
    }

    // Check second derivatives

    if (rbf_function->derivative_available(2))
    {
        for (int d = 0; d < dimension; ++d)
        {
            double dd_basis = rbf_function->dd_basis(group,
                                                     d,
                                                     shape,
                                                     r,
                                                     r0,
                                                     direction);
            
            if (!ce::approx(dd_basis, expected_double_grad[d + dimension * d], tol))
            {
                cout << "dd_basis in dimension ";
                cout << d;
                cout << " failed for ";
                cout << test_case;
                cout << endl;
                cout << "\texpected: ";
                cout << expected_double_grad[d + dimension * d];
                cout << "\tcalculated: ";
                cout << dd_basis;
                cout << endl;
                checksum += 1;
            }
        }

        double expected_laplacian = 0;
        
        for (int d = 0; d < dimension; ++d)
        {
            expected_laplacian += expected_double_grad[d + dimension * d];
        }

        double laplacian = rbf_function->laplacian(group,
                                                   shape,
                                                   r,
                                                   r0,
                                                   direction);
        
        if (!ce::approx(laplacian, expected_laplacian, tol))
        {
            cout << "laplacian failed for ";
            cout << test_case;
            cout << endl;
            cout << "\texpected: ";
            cout << expected_laplacian;
            cout << "\tcalculated: ";
            cout << laplacian;
            cout << endl;
            checksum += 1;
        }
    }
    
    return checksum;
}

int main()
{
    int checksum = 0;

    // Test 1
    
    {
        int const dimension = 2;
        
        vector<double> const r = {4,
                                  -3};
        vector<double> const r0 = {-2,
                                   7};
        vector<double> const direction = {1 / sqrt(3.),
                                          sqrt(2.) / sqrt(3.)};
        
        double const shape = 2.0;

        shared_ptr<RBF> rbf
            = make_shared<Multiquadric_RBF>();
        shared_ptr<Distance> distance
            = make_shared<Cartesian_Distance>(dimension);

        // Regular RBF
        
        {
            string test_case = "standard rbf";
            
            shared_ptr<RBF_Function> rbf_function
                = make_shared<RBF_Function>(rbf,
                                            distance);
            
            double const expected_basis = sqrt(545.);
            vector<double> const expected_grad = {24. / sqrt(545.),
                                                  -8. * sqrt(5. / 109.)};
            vector<double> const expected_double_grad = {1604. / (545. * sqrt(545.)),
                                                         192. / (109. * sqrt(545.)),
                                                         192. / (109. * sqrt(545.)),
                                                         116. / (109. * sqrt(545.))};
            
            checksum += test_rbf_function(rbf_function,
                                          test_case,
                                          dimension,
                                          shape,
                                          expected_basis,
                                          expected_grad,
                                          expected_double_grad,
                                          r,
                                          r0,
                                          direction);
        }

        // Derivative (directional) RBF
        
        {
            string test_case = "directional rbf";
            
            shared_ptr<RBF_Function> rbf_function
                = make_shared<Derivative_RBF_Function>(rbf,
                                                       distance);
            
            double const expected_basis = 8. * (3 - 5 * sqrt(2.)) / sqrt(1635.);
            vector<double> const expected_grad = {4. * (401. + 240. * sqrt(2.)) / (545. * sqrt(1635.)),
                                                  4. * (48. + 29. * sqrt(2)) / (109. * sqrt(1635.))};
            vector<double> const expected_double_grad(dimension * dimension, 0);
            
            checksum += test_rbf_function(rbf_function,
                                          test_case,
                                          dimension,
                                          shape,
                                          expected_basis,
                                          expected_grad,
                                          expected_double_grad,
                                          r,
                                          r0,
                                          direction);
        }
    }
    
    return checksum;
}
