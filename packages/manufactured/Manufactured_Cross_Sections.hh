#ifndef Manufactured_Cross_Sections_hh
#define Manufactured_Cross_Sections_hh

#include <memory>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;

/*
  Pure virtual class for a cross section in a manufactured solution
*/

class Manufactured_Cross_Sections
{
public:
    
    Manufactured_Cross_Sections(std::shared_ptr<Angular_Discretization> angular,
                                std::shared_ptr<Energy_Discretization> energy);

    // Get cross sections (sigma_f is combined with sigma_s)
    virtual void get_cross_sections(std::vector<double> const &position,
                                    std::vector<double> &sigma_t,
                                    std::vector<double> &sigma_s) const = 0;
    
protected:

    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
};

#endif
