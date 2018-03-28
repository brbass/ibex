#ifndef VERA_Heat_Data_hh
#define VERA_Heat_Data_hh

#include <memory>

#include "Heat_Transfer_Data.hh"

typedef std::function<double(std::vector<double> const &)> VERA_Temperature;
class VERA_Transport_Result;

class VERA_Heat_Data : public Heat_Transfer_Data
{
public:

    VERA_Heat_Data(std::shared_ptr<VERA_Transport_Result> result,
                   std::shared_ptr<VERA_Temperature> weighting_temperature);

    virtual double conduction(std::vector<double> const &position) const override;
    virtual double convection(std::vector<double> const &position) const override;
    virtual double source(std::vector<double> const &position) const override;
    virtual double temperature_inf(std::vector<double> const &position) const override;

private:

    std::shared_ptr<VERA_Transport_Result> result_;
    std::shared_ptr<VERA_Temperature> weighting_temperature_;
    // int number_of_materials_;
    double temperature_inf_;
    double convection_;
    // std::vector<double> conduction_;
    // std::vector<double> interfaces_;
};

#endif
