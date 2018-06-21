#include "VERA_Heat_Data.hh"

#include <cmath>

#include "VERA_Solid_Geometry.hh"
#include "VERA_Transport_Result.hh"

using namespace std;

VERA_Heat_Data::
VERA_Heat_Data(bool include_crack,
               int heat_dimension,
               shared_ptr<VERA_Transport_Result> result,
               shared_ptr<VERA_Temperature> weighting_temperature):
    Heat_Transfer_Data(),
    include_crack_(include_crack),
    heat_dimension_(heat_dimension),
    result_(result),
    weighting_temperature_(weighting_temperature)
{
    // number_of_materials_ = 3;
    temperature_inf_ = 600.0;
    convection_ = 3.0;
    // conduction_
    //     = {0.0377, 0.0245, 0.1649};
    if (include_crack)
    {
        interfaces_
            = {0.4135, 0.418, 0.475};
    }
    else
    {
        interfaces_
            = {0.4096, 0.418, 0.475};
    }
}

double VERA_Heat_Data::
conduction(vector<double> const &position) const
{
    // Get weighting temperature at position
    Assert(position.size() == heat_dimension_);
    double radius;
    vector<double> temp_position;
    switch (heat_dimension_)
    {
    case 1:
        temp_position = {position[0], 0};
        radius = position[0];
        break;
    case 2:
        temp_position = position;
        radius = sqrt(position[0] * position[0] + position[1] * position[1]);
        break;
    default:
        AssertMsg(false, "dimension not found");
    }
    double temperature = (*weighting_temperature_)(temp_position);

    // Get conduction
    if (include_crack_)
    {
        // Check for cracks
        Assert(false);
    }
    if (radius < interfaces_[0])
    {
        // Get fuel conduction
        double t = temperature / 1000;

        double k = 100. / (7.5408 + 17.692 * t + 3.6142 * t * t) + 6400. / pow(t, 2.5) * exp(-16.35 / t);
        return 0.01 * k;
    }
    else if (radius < interfaces_[1])
    {
        // Get gap conduction
        double k = 0.17632e-2 * pow(temperature, 0.77163);

        return 0.01 * k;
    }
    else
    {
        // Get clad conduction
        double t = temperature;
        double t2 = t * t;
        double t3 = t * t * t;
        double k = 7.51 + 2.09e-2 * t - 1.45e-5 * t2 + 7.69e-9 * t3;
        
        return 0.01 * k;
    }
}

double VERA_Heat_Data::
convection(vector<double> const &position) const
{
    return convection_;
}

double VERA_Heat_Data::
source(vector<double> const &position) const
{
    switch (heat_dimension_)
    {
    case 1:
        return result_->get_radial_fission_energy(position[0]);
    case 2:
        return result_->get_fission_energy(position);
    default:
        AssertMsg(false, "dimension not found");
        return -1;
    }
}

double VERA_Heat_Data::
temperature_inf(vector<double> const &position) const
{
    return temperature_inf_;
}
