#ifndef Manufactured_Heat_Data_hh
#define Manufactured_Heat_Data_hh

#include <memory>
#include <vector>

#include "Heat_Transfer_Data.hh"

/*
  Given conduction, convection and solution, returns the correct source and infinite medium temperature to give a certain solution
*/
class Manufactured_Heat_Data : public Heat_Transfer_Data
{
public:
    
    Manufactured_Heat_Data(int dimension,
                           std::vector<std::vector<double> > limits);
    
    // Functions to be implemented for each case
    virtual double conduction(std::vector<double> const &position) const override = 0;
    virtual std::vector<double> gradient_conduction(std::vector<double> const &position) const = 0;
    virtual double convection(std::vector<double> const &position) const override = 0;
    virtual double solution(std::vector<double> const &position) const = 0;
    virtual std::vector<double> gradient_solution(std::vector<double> const &position) const = 0;
    virtual double laplacian_solution(std::vector<double> const &position) const = 0;
    
    // Functions defined by the manufactured solution
    virtual double source(std::vector<double> const &position) const override;
    virtual double temperature_inf(std::vector<double> const &position) const override;
    virtual std::vector<double> normal(std::vector<double> const &position) const;
    
protected:
    
    int dimension_;
    std::vector<std::vector<double> > limits_;
};

#endif
