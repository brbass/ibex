#ifndef Material_Factory_hh
#define Material_Factory_hh

#include <memory>
#include <vector>

class Angular_Discretization;
class Energy_Discretization;
class Material;

/*
  Create a Material object standard data
*/
class Material_Factory
{
public:

    // Constructor
    Material_Factory(std::shared_ptr<Angular_Discretization> angular,
                     std::shared_ptr<Energy_Discretization> energy);
    
    // Return Material object
    std::shared_ptr<Material> get_standard_material(int index,
                                                    std::vector<double> const &sigma_t_data,
                                                    std::vector<double> const &sigma_s_data,
                                                    std::vector<double> const &nu_data,
                                                    std::vector<double> const &sigma_f_data,
                                                    std::vector<double> const &chi_data,
                                                    std::vector<double> const &internal_source_data) const;
    
private:
    
    std::shared_ptr<Angular_Discretization> angular_;
    std::shared_ptr<Energy_Discretization> energy_;
};

#endif
