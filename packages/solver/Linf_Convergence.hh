#ifndef Linf_Convergence_hh
#define Linf_Convergence_hh

#include "Convergence_Measure.hh"

class Linf_Convergence : public Convergence_Measure
{
public:

    Linf_Convergence(double tolerance = 1e-8);

    virtual void set_tolerance(double tolerance) override
    {
        tolerance_ = tolerance;
    }
    virtual double tolerance() const override
    {
        return tolerance_;
    }
    virtual double error(double val_new,
                         double val_old) const override;
    virtual double error(std::vector<double> const &val_new,
                         std::vector<double> const &val_old) const override;
    virtual bool check(double error_new,
                       double error_old) const override;

private:
    double tolerance_;
};

#endif
