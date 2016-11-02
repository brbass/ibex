#ifndef Transport_Discretization_hh
#define Transport_Discretization_hh

#include <memory>

class Angular_Discretization;
class Energy_Discretization;
class Spatial_Discretization;

using std::shared_ptr;

class Transport_Discretization
{
public:
    
    Transport_Discretization(shared_ptr<Spatial_Discretization> spatial,
                             shared_ptr<Angular_Discretization> angular,
                             shared_ptr<Energy_Discretization> energy);

    bool has_reflection() const
    {
        return has_reflection_;
    }
    int phi_size() const
    {
        return phi_size_;
    }
    int psi_size() const
    {
        return psi_size_;
    }
    int number_of_augments() const
    {
        return number_of_augments_;
    }
    
private:

    bool has_reflection_;
    int phi_size_;
    int psi_size_;
    int number_of_augments_;

    shared_ptr<Spatial_Discretization> spatial_;
    shared_ptr<Angular_Discretization> angular_;
    shared_ptr<Energy_Discretization> energy_;
};

#endif
