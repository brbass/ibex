#ifndef Convergence_Measure_hh
#define Convergence_Measure_hh

#include <vector>

class Convergence_Measure
{
public:
    Convergence_Measure();

    virtual double tolerance() const = 0;
    virtual double error(double val_new,
                         double val_old) const = 0;
    virtual double error(std::vector<double> const &val_new,
                         std::vector<double> const &val_old) const = 0;
    virtual bool check(double error_new,
                       double error_old) const = 0;
};

#endif
