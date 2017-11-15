#ifndef Quadratic_MLS_Normalization_hh
#define Quadratic_MLS_Normalization_hh

#include "Meshless_Normalization.hh"

#include <memory>

template<class T> class Dense_Solver;

class Quadratic_MLS_Normalization : public Meshless_Normalization
{
public:

    // Constructor
    Quadratic_MLS_Normalization(int dimension);

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
    void get_d_polynomial(int dim,
                          std::vector<double> const &position,
                          std::vector<double> &d_poly) const;
    void get_grad_polynomial(std::vector<double> const &position,
                             std::vector<std::vector<double> > &grad_poly) const;
    
    

private:

    // Data
    int dimension_;
    int number_of_polynomials_;
    std::shared_ptr<Dense_Solver<double> > solver_;
};


#endif
