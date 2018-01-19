#ifndef Manufactured_Sinusoidal_Cross_Sections_hh
#define Manufactured_Sinusoidal_Cross_Sections_hh

#include <memory>
#include <vector>

#include "Manufactured_Cross_Sections.hh"

class Angular_Discretization;
class Energy_Discretization;

/*
  Sinosoidal cross section of the form
  \Sigma = \Sigma_0 (1 + u \cos (2\pi f_x x) \cos (2\pi f_y y) \cos (2\pi f_z z) )
  f_d: frequency
  u: relative amplitude
*/

class Manufactured_Sinusoidal_Cross_Sections : public Manufactured_Cross_Sections
{
public:
    
    Manufactured_Sinusoidal_Cross_Sections(std::shared_ptr<Angular_Discretization> angular,
                                           std::shared_ptr<Energy_Discretization> energy,
                                           double const relative_amplitude,
                                           std::vector<double> const &frequency,
                                           std::vector<double> const &sigma_t,
                                           std::vector<double> const &sigma_s);

    // Get cross sections (sigma_f is combined with sigma_s)
    virtual void get_cross_sections(std::vector<double> const &position,
                                    std::vector<double> &sigma_t,
                                    std::vector<double> &sigma_s) const override;
    
private:
    
    double relative_amplitude_;
    std::vector<double> frequency_;
    std::vector<double> sigma_t_;
    std::vector<double> sigma_s_;
};

#endif
