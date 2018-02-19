#ifndef Strong_Spatial_Discretization_hh
#define Strong_Spatial_Discretization_hh

#include "Weak_Spatial_Discretization.hh"

class Strong_Spatial_Discretization : public Weak_Spatial_Discretization
{
public:
    
    // Constructor
    Strong_Spatial_Discretization(std::vector<std::shared_ptr<Basis_Function> > &bases,
                                  std::vector<std::shared_ptr<Weight_Function> > &weights,
                                  std::shared_ptr<Dimensional_Moments> dimensional_moments,
                                  std::shared_ptr<Weak_Spatial_Discretization_Options> options,
                                  std::shared_ptr<KD_Tree> kd_tree = std::shared_ptr<KD_Tree>());
private:
    void perform_basis_integration();
    void perform_point_integration();
};

#endif
