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

    cout << sum << endl;

    if (ce::approx(2 * (dimension - 1) * M_PI, sum, 1e-7))
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

int main()
{
    int checksum = 0;
    
    vector<int> ldfe_cases = {1, 2, 3};
    vector<int> ldfe_dims = {2, 3};
    int num_ldfe_cases = ldfe_cases.size();
    int num_ldfe_dims = ldfe_dims.size();
    
    for (int i = 0; i < num_ldfe_cases; ++i)
    {
        for (int j = 0; j < num_ldfe_dims; ++j)
        {
            checksum += test_ldfe(ldfe_dims[j],
                                  ldfe_cases[i]);
        }
    }

    return checksum;
}
