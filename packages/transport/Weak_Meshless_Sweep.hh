#ifndef Weak_Meshless_Sweep_hh
#define Weak_Meshless_Sweep_hh

#include "Meshless_Sweep.hh"

class Weak_Meshless_Sweep : public Meshless_Sweep
{
public:

    // Constructor
    Weak_Meshless_Sweep(Options options,
                        std::shared_ptr<Weak_Spatial_Discretization> spatial_discretization,
                        std::shared_ptr<Angular_Discretization> angular_discretization,
                        std::shared_ptr<Energy_Discretization> energy_discretization,
                        std::shared_ptr<Transport_Discretization> transport_discretization);

    virtual void check_class_invariants() const override;
    virtual std::string description() const override
    {
        return "Weak_Meshless_Sweep";
    }
    
protected:
    
    virtual void get_matrix_row(int i, // weight function index (row)
                                int o, // ordinate
                                int g, // group
                                std::vector<int> &indices, // global basis (column indices)
                                std::vector<double> &values) const override; // column values
    virtual void get_prec_matrix_row(int i, // weight function index (row)
                                     std::vector<int> &indices, // global basis (column indices)
                                     std::vector<double> &values) const override; // column values
    virtual void get_rhs(int i, // weight function index (row)
                         int o, // ordinate
                         int g, // group
                         std::vector<double> const &x, // angular flux w/ augments
                         double &value) const override; // rhs value
};

#endif
