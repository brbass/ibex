#ifndef Meshless_Normalization_hh
#define Meshless_Normalization_hh

#include <vector>

class Meshless_Normalization
{
public:
    
    Meshless_Normalization();

    // Should allow for base values and values to be same vector
    virtual void get_values(std::vector<double> const &position,
                            std::vector<std::vector<double> > const &center_positions,
                            std::vector<double> const &base_values,
                            std::vector<double> &values) const = 0;
    virtual void get_gradient_values(std::vector<double> const &position,
                                     std::vector<std::vector<double> > const &center_positions,
                                     std::vector<double> const &base_values,
                                     std::vector<std::vector<double> > const &grad_base_values,
                                     std::vector<double> &values,
                                     std::vector<std::vector<double> > &gradient_values) const = 0;
};

#endif
