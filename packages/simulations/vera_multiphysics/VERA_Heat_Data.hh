#ifndef VERA_Heat_Data_hh
#define VERA_Heat_Data_hh

#include "Heat_Transfer_Data.hh"

class VERA_Transport_Result;

class VERA_Heat_Data : public Heat_Transfer_Data
{
public:

    VERA_Heat_Data(std::shared_ptr<VERA_Transport_Result> result);

    virtual double conduction(std::vector<double> const &position) const override;
    virtual double convection(std::vector<double> const &position) const override;
    virtual double source(std::vector<double> const &position) const override;
    virtual double temperature_inf(std::vector<double> const &position) const override;

private:

    std::shared_ptr<VERA_Transport_Result> result_;
    int number_of_materials_;
    double temperature_inf_;
    double convection_;
    std::vector<double> conduction_;
    std::vector<double> interfaces_;
};

#endif
