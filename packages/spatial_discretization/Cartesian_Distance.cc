#include "Cartesian_Distance.hh"

#include <cmath>
#include <memory>

#include "Check.hh"

using std::make_shared;
using std::sqrt;
using std::vector;

Cartesian_Distance::
Cartesian_Distance(int dimension):
    dimension_(dimension)
{
}

double Cartesian_Distance::
distance(vector<double> const &r,
         vector<double> const &r0) const
{
    double val = 0;

    for (int i = 0; i < dimension_; ++i)
    {
        double const k = r[i] - r0[i];
        
        val += k * k;
    }
    
    return sqrt(val);
}

double Cartesian_Distance::
d_distance(int dim,
          vector<double> const &r,
          vector<double> const &r0) const
{
    double const dist = distance(r,
                                 r0);
    
    if (dist == 0.)
    {
        return 1;
    }
    
    double const distdim = r[dim] - r0[dim];
    
    return distdim / dist;
}

double Cartesian_Distance::
dd_distance(int dim,
            vector<double> const &r,
            vector<double> const &r0) const
{
    double const dist = distance(r,
                                 r0);
    
    if (dist == 0.)
    {
        return 0;
    }
    
    double const dist2 = dist * dist;
    double const dist3 = dist2 * dist;
    double const distdim = r[dim] - r0[dim];
    double const distdim2 = distdim * distdim;
    
    return (dist2 - distdim2) / dist3;
}

vector<double> Cartesian_Distance::
gradient_distance(vector<double> const &r,
                  vector<double> const &r0) const
{
    double const dist = distance(r,
                                 r0);
    
    vector<double> gradient(dimension_, 1);
    
    if (dist == 0.)
    {
        return gradient;
    }
    
    for (int d = 0; d < dimension_; ++d)
    {
        double const distdim = r[d] - r0[d];
        
        gradient[d] = distdim / dist;
    }
    
    return gradient;
}

vector<double> Cartesian_Distance::
double_gradient_distance(vector<double> const &r,
                         vector<double> const &r0) const
{
    vector<double> double_gradient(dimension_ * dimension_, 0);
    
    double const dist = distance(r,
                                 r0);
    
    if (dist == 0.)
    {
        return double_gradient;
    }
    
    vector<double> gradient = gradient_distance(r,
                                                r0);
    
    for (int d1 = 0; d1 < dimension_; ++d1)
    {
        for (int d2 = 0; d2 < dimension_; ++d2)
        {
            int k = d2 + dimension_ * d1;

            if (d1 == d2)
            {
                double_gradient[k] = (1 - gradient[d1] * gradient[d2]) / dist;
            }
            else
            {
                double_gradient[k] = -gradient[d1] * gradient[d2] / dist;
            }
        }
    }

    return double_gradient;
}

double Cartesian_Distance::
laplacian_distance(vector<double> const &r,
                   vector<double> const &r0) const
{
    double const dist = distance(r,
                                 r0);
    
    if (dist == 0.)
    {
        return 0;
    }
    
    return (dimension_ - 1) / dist;
}

void Cartesian_Distance::
check_class_invariants() const
{
    Assert(dimension_ > 0);
}
