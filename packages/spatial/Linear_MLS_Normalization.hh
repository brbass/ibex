#ifndef Linear_MLS_Normalization_hh
#define Linear_MLS_Normalization_hh

#include "Meshless_Normalization.hh"

class Linear_MLS_Normalization : public Meshless_Normalization
{
public:

    // Constructor
    Linear_MLS_Normalization(int dimension);

    // Methods from Meshless_Normalization
    virtual void get_values(std::vector<double> const &position,
                            std::vector<std::vector<double> > const &center_positions,
                            std::vector<double> const &base_values,
                            std::vector<double> &values) const override;
    virtual void get_gradient_values(std::vector<double> const &position,
                                     std::vector<std::vector<double> > const &center_positions,
                                     std::vector<double> const &base_values,
                                     std::vector<std::vector<double> > const &grad_base_values,
                                     std::vector<double> &values,
                                     std::vector<std::vector<double> > &gradient_values) const override;
    
    // Helper methods
    void get_polynomial(std::vector<double> const &position,
                        std::vector<double> &poly) const;
    void get_grad_polynomial(std::vector<double> const &position,
                             std::vector<std::vector<double> > &grad_poly) const;
    
    

private:

    // Data
    int dimension_;
    int number_of_polynomials_;
};


#endif
