#include <cmath>
#include <iostream>

#include "Check_Equality.hh"
#include "LDFE_Quadrature.hh"

namespace ce = Check_Equality;

using namespace std;

int test_ldfe(int dimension,
              int rule)
{
    int number_of_moments = 1;

    LDFE_Quadrature quad(dimension,
                         number_of_moments,
                         rule);

    int number_of_ordinates = quad.number_of_ordinates();
    vector<double> const weights = quad.weights();
    
    double sum = 0;

    for (int i = 0; i < number_of_ordinates; ++i)
    {
        sum += weights[i];
    }

    if (ce::approx(2 * (dimension - 1) * M_PI, sum, 1e-7))
    {
        return 0;
    }
    else
    {
        cout << sum << endl;
        return 1;
    }
}

int main()
{
    int checksum = 0;
    
    vector<int> dimensions = {2, 3};
    vector<int> rules = {1, 2, 3, 4, 5, 6};

    for (int dimension : dimensions)
    {
        for (int rule : rules)
        {
            checksum += test_ldfe(dimension,
                                  rule);
        }
    }
    
    return checksum;
}
