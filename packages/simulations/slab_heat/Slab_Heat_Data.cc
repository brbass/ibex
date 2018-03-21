#include "Slab_Heat_Data.hh"

#include <cmath>
#include <limits>

#include "Check.hh"
#include "Heat_Transfer_Integration.hh"

using namespace std;

Slab_Heat_Data::
Slab_Heat_Data(shared_ptr<Heat_Transfer_Integration_Options> int_options,
               std::vector<double> k,
               std::vector<double> q,
               std::vector<double> h,
               std::vector<double> tinf,
               std::vector<double> xlim):
    int_options_(int_options),
    k_(k),
    q_(q),
    h_(h),
    tinf_(tinf),
    xlim_(xlim)
{
    boundary_tol_ = 100 * numeric_limits<double>::epsilon();
    AssertMsg(abs(8.*k_[0]*k_[1]*(h_[1]*k_[0] + h_[0]*k_[1])*pow(q_[2],2)) > boundary_tol_,
              "parameters incompatible with benchmark solution");
}

double Slab_Heat_Data::
conduction(vector<double> const &position) const
{
    double x = position[0];
    return x < 0 ? k_[0] : k_[1];
}

double Slab_Heat_Data::
convection(vector<double> const &position) const
{
    if (abs(position[0] - xlim_[0]) < boundary_tol_)
    {
        return h_[0];
    }
    else if (abs(position[0] - xlim_[1]) < boundary_tol_)
    {
        return h_[1];
    }
    else
    {
        return 0;
    }
}

double Slab_Heat_Data::
source(vector<double> const &position) const
{
    double x = position[0];
    return q_[0] + pow(sin(q_[1] * x), 2);
}

double Slab_Heat_Data::
temperature_inf(vector<double> const &position) const
{
    double x = position[0];
    return x < 0 ? tinf_[0] : tinf_[1];
}

double Slab_Heat_Data::
get_solution(vector<double> const &position) const
{
    double x = position[0];
    if (x < 0)
    {
        return (-(cos(2*x*q_[2])*(h_[1]*k_[0] + h_[0]*k_[1])*q_[1]) + cos(2*q_[2]*xlim_[0])*(h_[1]*k_[0] + h_[0]*k_[1])*q_[1] + 2*q_[2]*(-(k_[0]*q_[1]*sin(2*q_[2]*xlim_[1])*(k_[0] + h_[0]*(x - xlim_[0]))) + k_[0]*q_[1]*sin(2*q_[2]*xlim_[0])*(k_[1] + h_[1]*(-x + xlim_[0])) + q_[2]*(-(pow(x,2)*(h_[1]*k_[0] + h_[0]*k_[1])*(2*q_[0] + q_[1])) + h_[0]*(k_[1]*(2*q_[0] + q_[1])*pow(xlim_[0],2) + 4*k_[0]*(k_[1]*tinf_[0] + h_[1]*(tinf_[0] - tinf_[1])*xlim_[0]) - 2*k_[0]*(2*q_[0] + q_[1])*xlim_[0]*xlim_[1]) + k_[0]*(h_[1]*(4*k_[0]*tinf_[1] - (2*q_[0] + q_[1])*pow(xlim_[0],2)) - 2*(2*q_[0] + q_[1])*(k_[1]*xlim_[0] - k_[0]*xlim_[1])) + 2*x*k_[0]*(h_[1]*(2*q_[0] + q_[1])*xlim_[0] + h_[0]*(2*h_[1]*(-tinf_[0] + tinf_[1]) + (2*q_[0] + q_[1])*xlim_[1])))))/(8.*k_[0]*(h_[1]*k_[0] + h_[0]*k_[1])*pow(q_[2],2));
    }
    else
    {
        return (-(cos(2*x*q_[2])*k_[0]*(h_[1]*k_[0] + h_[0]*k_[1])*q_[1]) + cos(2*q_[2]*xlim_[0])*k_[1]*(h_[1]*k_[0] + h_[0]*k_[1])*q_[1] - (h_[1]*k_[0] + h_[0]*k_[1])*((-k_[0] + k_[1])*q_[1] + 2*pow(x,2)*k_[0]*(2*q_[0] + q_[1])*pow(q_[2],2)) + 2*k_[1]*q_[2]*(-(k_[0]*q_[1]*sin(2*q_[2]*xlim_[1])*(k_[0] + h_[0]*(x - xlim_[0]))) + k_[0]*q_[1]*sin(2*q_[2]*xlim_[0])*(k_[1] + h_[1]*(-x + xlim_[0])) + q_[2]*(h_[0]*(k_[1]*(2*q_[0] + q_[1])*pow(xlim_[0],2) + 4*k_[0]*(k_[1]*tinf_[0] + h_[1]*(tinf_[0] - tinf_[1])*xlim_[0]) - 2*k_[0]*(2*q_[0] + q_[1])*xlim_[0]*xlim_[1]) + k_[0]*(h_[1]*(4*k_[0]*tinf_[1] - (2*q_[0] + q_[1])*pow(xlim_[0],2)) - 2*(2*q_[0] + q_[1])*(k_[1]*xlim_[0] - k_[0]*xlim_[1])) + 2*x*k_[0]*(h_[1]*(2*q_[0] + q_[1])*xlim_[0] + h_[0]*(2*h_[1]*(-tinf_[0] + tinf_[1]) + (2*q_[0] + q_[1])*xlim_[1])))))/(8.*k_[0]*k_[1]*(h_[1]*k_[0] + h_[0]*k_[1])*pow(q_[2],2));
    }
}
