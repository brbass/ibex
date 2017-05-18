#ifndef Linear_MLS_External_hh
#define Linear_MLS_External_hh

#include <vector>

class Linear_MLS_External
{
public:
    
    Linear_MLS_External(int dimension);

    void get_polynomial(std::vector<double> const &position,
                        std::vector<double> &poly) const;
    void get_grad_polynomial(std::vector<double> const &position,
                             std::vector<std::vector<double> > &grad_poly) const;
    
    
    void get_values(std::vector<double> const &position,
                    std::vector<std::vector<double> > const &center_positions,
                    std::vector<double> const &base_values,
                    std::vector<double> &values) const;
    void get_gradient_values(std::vector<double> const &position,
                             std::vector<std::vector<double> > const &center_positions,
                             std::vector<double> const &base_values,
                             std::vector<std::vector<double> > const &grad_base_values,
                             std::vector<double> &values,
                             std::vector<std::vector<double> > &gradient_values) const;
    
private:

    int dimension_;
    int number_of_polynomials_;
};


#endif
