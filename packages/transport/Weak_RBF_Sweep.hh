#ifndef Weak_RBF_Sweep_hh
#define Weak_RBF_Sweep_hh

#include "Sweep_Operator.hh"
#include "Weak_Spatial_Discretization.hh"

class Weak_RBF_Sweep : public Sweep_Operator
{
public:
    
    struct Options
    {
        enum class Solver
        {
            AMESOS,
            AZTEC,
            EIGEN
        };

        bool use_supg;
        Solver solver = Solver::AMESOS;
    };

    // Constructor
    Weak_RBF_Sweep(Options options,
                   std::shared_ptr<Weak_Spatial_Discretization> spatial_discretization,
                   std::shared_ptr<Angular_Discretization> angular_discretization,
                   std::shared_ptr<Energy_Discretization> energy_discretization,
                   std::shared_ptr<Transport_Discretization> transport_discretization);

    // Sweep_Operator functions
    virtual std::shared_ptr<Spatial_Discretization> spatial_discretization() const override
    {
        return spatial_discretization;
    }
    virtual std::shared_ptr<Angular_Discretization> angular_discretization() const override
    {
        return angular_discretization;
    }
    virtual std::shared_ptr<Energy_Discretization> energy_discretization() const override
    {
        return energy_discretization;
    }
    virtual void output(XML_Node output_node) const override;
    virtual void check_class_invariants() const override;
    
private:

    // Vector_Operator function
    virtual void apply(std::vector<double> &x) const override;

    // Weak_RBF_Sweep functions
    void get_matrix_row(int i, // row
                        int o, // ordinate
                        int g, // group
                        vector<int> &indices,
                        vector<double> &values) const;
    
    // Data
    std::shared_ptr<Weak_Spatial_Discretization> spatial_discretization_;
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;
    std::shared_ptr<Transport_Discretization> transport_discretization_;
};

#endif
