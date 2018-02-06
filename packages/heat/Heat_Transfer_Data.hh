#ifndef Heat_Transfer_Data_hh
#define Heat_Transfer_Data_hh

#include <vector>

class Heat_Transfer_Data
{
public:

    Heat_Transfer_Data();

    virtual double conduction(std::vector<double> const &position) const = 0;
    virtual double convection(std::vector<double> const &position) const = 0;
    virtual double source(std::vector<double> const &position) const = 0;
    virtual double temperature_inf(std::vector<double> const &position) const = 0;
};

#endif
