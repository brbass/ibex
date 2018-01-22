#ifndef Manufactured_Slab_Cross_Sections_hh
#define Manufactured_Slab_Cross_Sections_hh

#include <memory>
#include <vector>

#include "Manufactured_Cross_Sections.hh"

class Angular_Discretization;
class Energy_Discretization;

/*
  Pure virtual class for a cross section in a manufactured solution
*/

class Manufactured_Slab_Cross_Sections : public Manufactured_Cross_Sections
{
public:
    
    Manufactured_Slab_Cross_Sections(std::shared_ptr<Angular_Discretization> angular,
                                     std::shared_ptr<Energy_Discretization> energy,
                                     std::vector<double> const &interface_positions,
                                     std::vector<std::vector<double> > const &sigma_t,
                                     std::vector<std::vector<double> > const &sigma_s);

    // Get cross sections (sigma_f is combined with sigma_s)
    virtual void get_cross_sections(std::vector<double> const &position,
                                    std::vector<double> &sigma_t,
                                    std::vector<double> &sigma_s) const override;
    
private:

    int number_of_regions_;
    std::vector<double> interface_positions_;
    std::vector<std::vector<double> > sigma_t_;
    std::vector<std::vector<double> > sigma_s_;
};

#endif
