#include "VERA_Heat_Data.hh"

#include "VERA_Transport_Result.hh"

using namespace std;

VERA_Heat_Data::
VERA_Heat_Data(shared_ptr<VERA_Transport_Result> result):
    Heat_Transfer_Data(),
    result_(result)
{
    number_of_materials_ = 3;
    temperature_inf_ = 600.0;
    convection_ = 3.0;
    conduction_
        = {0.0377, 0.0245, 0.1649};
    interfaces_
        = {0.4096, 0.418, 0.475};
}

double VERA_Heat_Data::
conduction(vector<double> const &position) const
{
    for (int i = 0; i < number_of_materials_; ++i)
    {
        if (position[0] < interfaces_[i])
        {
            return conduction_[i];
        }
    }

    return conduction_[number_of_materials_ - 1];
}

double VERA_Heat_Data::
convection(vector<double> const &position) const
{
    return convection_;
}

double VERA_Heat_Data::
source(vector<double> const &position) const
{
    return result_->get_radial_fission_energy(position[0]);
}

double VERA_Heat_Data::
temperature_inf(vector<double> const &position) const
{
    return temperature_inf_;
}
