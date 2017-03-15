#ifndef Multiplicative_Operator_hh
#define Multiplicative_Operator_hh

#include <vector>

#include "Vector_Operator.hh"

/*
  Return the given vector
*/
class Multiplicative_Operator : public Vector_Operator
{
public:

    // Constructor
    Multiplicative_Operator(int size,
                            double scalar);

    void set_scalar(double scalar)
    {
        scalar_ = scalar;
    }
    
    virtual void check_class_invariants() const override;

private:

    double scalar_;
    
    virtual void apply(std::vector<double> &x) const override;
};

#endif
