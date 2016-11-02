#ifndef RBF_Collocation_Sweep_hh
#define RBF_Collocation_Sweep_hh

#include <memory>

#include "RBF_Discretization.hh"
#include "Sweep_Operator.hh"

class Amesos_BaseSolver;
class AztecOO;
class Epetra_CrsMatrix;
class Epetra_Comm;
class Epetra_LinearProblem;
class Epetra_Map;
class Epetra_SerialDenseMatrix;
class Epetra_SerialDenseSolver;
class Epetra_Vector;

using std::shared_ptr;

class RBF_Collocation_Sweep : public Sweep_Operator
{
public:

    enum class Matrix_Solver
    {
        AMESOS,
        AZTEC
    };
    
    enum class Solution_Variable
    {
        COEFFICIENT,
        PSI,
    };

    enum class Condition_Calculation
    {
        NONE,
        CHEAP,
        EXPENSIVE,
        AZTEC
    };
    
    RBF_Collocation_Sweep(Solution_Variable solution_variable,
                          Matrix_Solver matrix_solver,
                          Condition_Calculation condition_calculation,
                          shared_ptr<RBF_Discretization> rbf_discretization,
                          shared_ptr<Angular_Discretization> angular_discretization,
                          shared_ptr<Energy_Discretization> energy_discretization,
                          shared_ptr<Transport_Discretization> transport_discretization);

    virtual shared_ptr<RBF_Discretization> rbf_discretization() const
    {
        return rbf_discretization_;
    }
    virtual shared_ptr<Spatial_Discretization> spatial_discretization() const override
    {
        return rbf_discretization_;
    }
    virtual shared_ptr<Angular_Discretization> angular_discretization() const override
    {
        return angular_discretization_;
    }
    virtual shared_ptr<Energy_Discretization> energy_discretization() const override
    {
        return energy_discretization_;
    }

    virtual void output(pugi::xml_node output_node) const override;
    virtual void check_class_invariants() const override;

private:

    virtual void apply(vector<double> &x) const override;
    
    void calculate_condition_numbers();
    void initialize_trilinos();
    void set_point(int i,
                   int o,
                   int g);
    void set_boundary_point(int i,
                            int o,
                            int g);
    void set_internal_point(int i,
                            int o,
                            int g);
    void set_boundary_rhs(int b,
                          int i,
                          int o,
                          int g,
                          vector<double> const &x) const;
    void set_internal_rhs(int i,
                          int o,
                          int g,
                          vector<double> const &x) const;

    void convert_to_psi(int i,
                        int g,
                        int o,
                        vector<double> &data) const;
    
    void update_local_matrix(int i,
                             int g,
                             int o);
    void check_local_matrix(int i,
                            int g,
                            int o);
    
    Solution_Variable solution_variable_;
    Condition_Calculation condition_calculation_;
    Matrix_Solver matrix_solver_;
    shared_ptr<RBF_Discretization> rbf_discretization_;
    shared_ptr<Angular_Discretization> angular_discretization_;
    shared_ptr<Energy_Discretization> energy_discretization_;

    double reflection_tolerance_;

    shared_ptr<Epetra_Comm> comm_;
    shared_ptr<Epetra_Map> map_;
    shared_ptr<Epetra_Vector> lhs_;
    shared_ptr<Epetra_Vector> rhs_;
    vector<shared_ptr<Epetra_CrsMatrix> > mat_;
    vector<shared_ptr<Epetra_LinearProblem> > problem_;

    int max_iterations_;
    double tolerance_;
    mutable vector<int> num_calls_;
    mutable vector<int> num_iterations_;
    vector<shared_ptr<AztecOO> > aztec_solver_;
    vector<shared_ptr<Amesos_BaseSolver*> > amesos_solver_;
    
    vector<double> condition_numbers_;
    
    shared_ptr<Epetra_SerialDenseMatrix> local_mat_;
    shared_ptr<Epetra_SerialDenseSolver> local_solver_;
    mutable vector<int> num_averages_;
    mutable vector<double> average_local_condition_numbers_;
};

#endif
