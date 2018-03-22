#include "Manufactured_Sin_Data.hh"

#include <cmath>
#include <limits>

#include "Check.hh"

using namespace std;

Manufactured_Sin_Data::
Manufactured_Sin_Data(int dimension,
                       vector<vector<double> > limits,
                       vector<double> sol_coeff,
                       vector<double> cond_coeff,
                       vector<double> conv_coeff):
    Manufactured_Heat_Data(dimension,
                           limits),
    sol_coeff_(sol_coeff),
    cond_coeff_(cond_coeff),
    conv_coeff_(conv_coeff)
{
    Assert(limits_.size() == dimension_);
    Assert(sol_coeff_.size() == 2);
    Assert(cond_coeff.size() == 2);
    Assert(conv_coeff.size() == 3);
}

double Manufactured_Sin_Data::
conduction(std::vector<double> const &position) const
{
    double k1 = cond_coeff_[0];
    double k2 = cond_coeff_[1];
    double dist2 = 0;
    for (int d = 0; d < dimension_; ++d)
    {
        dist2 += position[d] * position[d];
    }
    
    return k1 + k2 * dist2;
}

std::vector<double> Manufactured_Sin_Data::
gradient_conduction(std::vector<double> const &position) const
{
    double k2 = cond_coeff_[1];
    vector<double> val(dimension_);
    for (int d = 0; d < dimension_; ++d)
    {
        val[d] = 2 * k2 * position[d];
    }

    return val;
}

double Manufactured_Sin_Data::
convection(std::vector<double> const &position) const
{
    double h1 = conv_coeff_[0];
    double h2 = conv_coeff_[1];
    double h3 = conv_coeff_[2];

    double posprod = 1;
    for (int d = 0; d < dimension_; ++d)
    {
        posprod *= position[d];
    }
    
    return h1 + h2 * sin(h3 * posprod);
}

double Manufactured_Sin_Data::
solution(std::vector<double> const &position) const
{
    double t1 = sol_coeff_[0];
    double t2 = sol_coeff_[1];

    double dist2 = 0;
    for (int d = 0; d < dimension_; ++d)
    {
        dist2 += position[d] * position[d];
    }

    return t1 + sin(t2 * dist2);
}

std::vector<double> Manufactured_Sin_Data::
gradient_solution(std::vector<double> const &position) const
{
    double t2 = sol_coeff_[1];
    
    double dist2 = 0;
    for (int d = 0; d < dimension_; ++d)
    {
        dist2 += position[d] * position[d];
    }

    double cosval = cos(t2 * dist2);
    vector<double> val(dimension_);
    for (int d = 0; d < dimension_; ++d)
    {
        val[d] = 2 * t2 * position[d] * cosval;
    }

    return val;
}

double Manufactured_Sin_Data::
laplacian_solution(std::vector<double> const &position) const
{
    double t2 = sol_coeff_[1];
    
    double dist2 = 0;
    for (int d = 0; d < dimension_; ++d)
    {
        dist2 += position[d] * position[d];
    }

    double cosval = cos(t2 * dist2);
    double sinval = sin(t2 * dist2);
    double val = 0;
    for (int d = 0; d < dimension_; ++d)
    {
        val += 2 * t2 * (cosval - 2 * t2 * position[d] * position[d] * sinval);
    }

    return val;
}

