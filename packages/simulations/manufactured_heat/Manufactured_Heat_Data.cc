#include "Manufactured_Heat_Data.hh"

#include <cmath>
#include <limits>

#include "Check.hh"

using namespace std;

Manufactured_Heat_Data::
Manufactured_Heat_Data(int dimension,
                       std::vector<std::vector<double> > limits):
    dimension_(dimension),
    limits_(limits)
{
    
}

double Manufactured_Heat_Data::
source(vector<double> const &position) const
{
    double k = conduction(position);
    vector<double> k_grad = gradient_conduction(position);
    vector<double> t_grad = gradient_solution(position);
    double t_lap = laplacian_solution(position);

    double val = -k * t_lap;
    for (int d = 0; d < dimension_; ++d)
    {
        val -= k_grad[d] * t_grad[d];
    }
    
    return val;
}

double Manufactured_Heat_Data::
temperature_inf(vector<double> const &position) const
{
    double h = convection(position);
    vector<double> n = normal(position);
    double k = conduction(position);
    double t = solution(position);
    vector<double> t_grad = gradient_solution(position);

    double val = 0;
    for (int d = 0; d < dimension_; ++d)
    {
        val += k * n[d] * t_grad[d];
    }
    val /= h;
    val += t;

    return val;
}

vector<double> Manufactured_Heat_Data::
normal(vector<double> const &position) const
{
    double const tolerance = 1000 * numeric_limits<double>::epsilon();
    
    vector<double> normal(dimension_, 0);
    for (int d = 0; d < dimension_; ++d)
    {
        for (int b = 0; b < 2; ++b)
        {
            if (abs(position[d] - limits_[d][b]) < tolerance)
            {
                normal[d] = b == 0 ? -1 : 1;
                
                return normal;
            }
        }
    }

    AssertMsg(false, "no boundary found at specified location");
    
    return normal;
}
