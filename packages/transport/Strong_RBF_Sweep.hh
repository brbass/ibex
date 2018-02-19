#ifndef Strong_RBF_Sweep_hh
#define Strong_RBF_Sweep_hh

#include "Weak_RBF_Sweep.hh"

class Strong_RBF_Sweep : public Weak_RBF_Sweep
{
public:
    
    Strong_RBF_Sweep(Options options,
                     std::shared_ptr<Weak_Spatial_Discretization> spatial_discretization,
                     std::shared_ptr<Angular_Discretization> angular_discretization,
                     std::shared_ptr<Energy_Discretization> energy_discretization,
                     std::shared_ptr<Transport_Discretization> transport_discretization);

protected:

    virtual void get_matrix_row(int i, // weight function index (row)
                                int o, // ordinate
                                int g, // group
                                std::vector<int> &indices, // global basis (column indices)
                                std::vector<double> &values) const override; // column values
    virtual void get_rhs(int i, // weight function index (row)
                         int o, // ordinate
                         int g, // group
                         std::vector<double> const &x, // angular flux w/ augments
                         double &value) const override; // rhs value
};

#endif
