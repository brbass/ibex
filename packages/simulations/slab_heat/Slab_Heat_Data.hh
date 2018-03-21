#ifndef Slab_Heat_Data_hh
#define Slab_Heat_Data_hh

#include <memory>
#include <vector>

#include "Heat_Transfer_Data.hh"

class Heat_Transfer_Integration_Options;

/*
  Represents a problem with:
  k(x,y,z) = {k[0], x < 0; k[0], otherwise}
  q(x,y,z) = q[0] + q[1] * sin^2(q[2] * x)
  h(x,y,z) = {h[0], x==xlim[0]; h[1], x==xlim[1]; 0, otherwise}
  tinf(x,y,z) = {tinf[0], x==xlim[0]; tinf[1], x==xlim[1]; 0, otherwise}
*/
class Slab_Heat_Data : public Heat_Transfer_Data
{
public:
    
    Slab_Heat_Data(std::shared_ptr<Heat_Transfer_Integration_Options> int_options,
                   std::vector<double> k,
                   std::vector<double> q,
                   std::vector<double> h,
                   std::vector<double> tinf,
                   std::vector<double> xlim);
    
    virtual double conduction(std::vector<double> const &position) const override;
    virtual double convection(std::vector<double> const &position) const override;
    virtual double source(std::vector<double> const &position) const override;
    virtual double temperature_inf(std::vector<double> const &position) const override;

    double get_solution(std::vector<double> const &position) const;
    
private:

    std::shared_ptr<Heat_Transfer_Integration_Options> int_options_;
    double boundary_tol_;
    std::vector<double> k_;
    std::vector<double> q_;
    std::vector<double> h_;
    std::vector<double> tinf_;
    std::vector<double> xlim_;
};

#endif
