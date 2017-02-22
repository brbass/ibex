#include <cmath>
#include <iostream>
#include <limits>

#include "Check_Equality.hh"
#include "Gauss_Legendre_Quadrature.hh"

namespace ce = Check_Equality;

using namespace std;

// x^(2o-1) + x^(2o-2)
double test_function1(int number_of_ordinates,
                      double x)
{
    int n = number_of_ordinates * 2 - 1;
    
    return pow(x, n) + pow(x, n - 1);
}
double integral_test_function1(int number_of_ordinates)
{
    int n = number_of_ordinates * 2 - 1;

    return (1 + pow(-1, n + 1) + 2 * n) / (n * (1 + n));
}

// \sum_{i=0}^{2o-1} n x^n
double test_function2(int number_of_ordinates,
                      double x)
{
    int n = number_of_ordinates * 2 - 1;

    return (x + pow(x, n + 1) * (n * x - n - 1)) / pow(1 - x, 2);
}
double integral_test_function2(int number_of_ordinates)
{
    int n = number_of_ordinates * 2 - 1;
    
    double result = 0;

    for (int i = n; i >= 0; --i)
    {
        result += (1 + pow(-1, i)) * i / (1 + i);
    }
    
    return result;
}

int test_gauss_legendre(int number_of_ordinates)
{
    int checksum = 0;
    
    int dimension = 1;
    int number_of_moments = 1;

    double tolerance = 1000 * numeric_limits<double>::epsilon();
    
    Gauss_Legendre_Quadrature quad(dimension,
                                   number_of_moments,
                                   number_of_ordinates);

    vector<double> const ordinates = quad.ordinates();
    vector<double> const weights = quad.weights();

    // Check that weights sum to 2
    
    double sum = 0;
    
    for (int i = 0; i < number_of_ordinates; ++i)
    {
        sum += weights[i];
    }
    
    if (!ce::approx(2.0, sum, tolerance))
    {
        cout << "gauss_legendre weights failed for order ";
        cout << number_of_ordinates;
        cout << endl;
        checksum += 1;
    }
    
    // Check that ordinates correctly integrate polynomial

    int number_of_functions = 2;
    for (int f = 0; f < number_of_functions; ++f)
    {
        double actual_integral;
        double integral = 0;
        double local_tolerance;
        switch(f)
        {
        case 0:
            actual_integral = integral_test_function1(number_of_ordinates);
            local_tolerance = tolerance;
            for (int i = 0; i < number_of_ordinates; ++i)
            {
                integral += weights[i] * test_function1(number_of_ordinates,
                                                        ordinates[i]);
            }
            break;
        case 1:
            actual_integral = integral_test_function2(number_of_ordinates);
            local_tolerance = actual_integral * number_of_ordinates * tolerance;
            for (int i = 0; i < number_of_ordinates; ++i)
            {
                integral += weights[i] * test_function2(number_of_ordinates,
                                                        ordinates[i]);
            }
            break;
        }
        
        if (!ce::approx(actual_integral,
                        integral,
                        local_tolerance))
        {
            cout << "gauss_legendre did not integrate correctly for order ";
            cout << number_of_ordinates;
            cout << " for function ";
            cout << f;
            cout << endl;
            cout << "\texpected:   ";
            cout << actual_integral;
            cout << "\tcalculated: ";
            cout << integral;
            cout << endl;
        }
    }
    return checksum;
    
}

int main()
{
    int checksum = 0;

    int num_gauss_cases = 100;
    
    for (int i = 0; i < num_gauss_cases; ++i)
    {
        int number_of_ordinates = 2 * (i + 1);
        checksum += test_gauss_legendre(number_of_ordinates);
    }
    
    return checksum;
}
