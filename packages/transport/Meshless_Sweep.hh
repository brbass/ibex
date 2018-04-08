#ifndef Meshless_Sweep_hh
#define Meshless_Sweep_hh

#include "Sweep_Operator.hh"
#include "Weak_Spatial_Discretization.hh"

class Amesos_BaseSolver;
class AztecOO;
template<class T1, class T2> class Conversion;
class Epetra_CrsMatrix;
class Epetra_Comm;
class Epetra_LinearProblem;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Operator;
class Epetra_Vector;
class Ifpack_Preconditioner;

namespace Belos
{
    class EpetraPrecOp;
    template<class Scalar, class MV, class OP> class LinearProblem;
    template<class Scalar, class MV, class OP> class PseudoBlockGmresSolMgr;
}

typedef Belos::EpetraPrecOp BelosPreconditioner;
typedef Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator> BelosLinearProblem;
typedef Belos::PseudoBlockGmresSolMgr<double, Epetra_MultiVector, Epetra_Operator> BelosSolver;

class Meshless_Sweep : public Sweep_Operator
{
public:

    // Sweep options
    struct Options
    {
        // Solver type
        enum class Solver
        {
            AMESOS,
            AMESOS_PARALLEL,
            AZTEC,
            AZTEC_IFPACK,
            BELOS,
            BELOS_IFPACK,
            BELOS_IFPACK_RIGHT,
            BELOS_IFPACK_RIGHT2
        };
        std::shared_ptr<Conversion<Solver, std::string> > solver_conversion() const;
        
        Solver solver = Solver::BELOS_IFPACK;

        // List of possible solver options
        bool quit_if_diverged = true;
        bool use_preconditioner = true;
        bool print = false;
        int max_iterations = 1000;
        int kspace = 20;
        int max_restarts = 50;
        double level_of_fill = 1.0;
        double tolerance = 1e-8;
        double drop_tolerance = 1e-12;

        // Options specific to right preconditioners
        bool weighted_preconditioner = false; 
        bool force_left = false;
    };

    // Constructor
    Meshless_Sweep(Options options,
                   std::shared_ptr<Weak_Spatial_Discretization> spatial_discretization,
                   std::shared_ptr<Angular_Discretization> angular_discretization,
                   std::shared_ptr<Energy_Discretization> energy_discretization,
                   std::shared_ptr<Transport_Discretization> transport_discretization);

    // Sweep_Operator functions
    virtual std::shared_ptr<Spatial_Discretization> spatial_discretization() const override
    {
        return spatial_discretization_;
    }
    virtual std::shared_ptr<Angular_Discretization> angular_discretization() const override
    {
        return angular_discretization_;
    }
    virtual std::shared_ptr<Energy_Discretization> energy_discretization() const override
    {
        return energy_discretization_;
    }
    virtual void output(XML_Node output_node) const override;
    virtual void check_class_invariants() const override = 0;
    virtual std::string description() const override = 0;
    
    // Save matrix to specified XML output file
    void save_matrix_as_xml(int o,
                            int g,
                            XML_Node output_node) const;
    
protected:

    // Vector_Operator function
    virtual void apply(std::vector<double> &x) const override;

    // Meshless_Sweep functions
    virtual void initialize_solver();
    virtual void update_augments(std::vector<double> &x) const;
    virtual void get_matrix_row(int i, // weight function index (row)
                                int o, // ordinate
                                int g, // group
                                std::vector<int> &indices, // global basis (column indices)
                                std::vector<double> &values) const = 0; // column values
    virtual void get_prec_matrix_row(int i, // weight function index (row)
                                     std::vector<int> &indices, // global basis (column indices)
                                     std::vector<double> &values) const = 0; // column values
    virtual void get_rhs(int i, // weight function index (row)
                         int o, // ordinate
                         int g, // group
                         std::vector<double> const &x, // angular flux w/ augments
                         double &value) const = 0; // rhs value
    
    // Generalized solver
    class Sweep_Solver
    {
    public:
        // Constructor
        Sweep_Solver(Meshless_Sweep const &wrs);

        // Solve problem
        virtual void solve(std::vector<double> &x) const = 0;

    protected:
        // Data
        Meshless_Sweep const &wrs_;
    };
    
    // Generalized trilinos solver
    class Trilinos_Solver : public Sweep_Solver
    {
    public:
        // Constructor
        Trilinos_Solver(Meshless_Sweep const &wrs);

        // Solve problem
        virtual void solve(std::vector<double> &x) const = 0;

    protected:
        
        // Return transport matrix for o and g
        std::shared_ptr<Epetra_CrsMatrix> get_matrix(int o,
                                                     int g,
                                                     std::shared_ptr<Epetra_Map> map) const;

        // Get preconditioner matrix that is independent of o and g
        std::shared_ptr<Epetra_CrsMatrix> get_prec_matrix(std::shared_ptr<Epetra_Map> map) const;
        
        // Set rhs_ to current source
        void set_rhs(int o,
                     int g,
                     std::shared_ptr<Epetra_Vector> &rhs,
                     std::vector<double> const &x) const;
        
        // Check Aztec solver message
        void check_aztec_convergence(std::shared_ptr<AztecOO> const solver) const;
    };
    
    // Amesos solver
    // Stores LU decompositions of all matrices
    // Only serial
    class Amesos_Solver : public Trilinos_Solver
    {
    public:
        // Constructor
        Amesos_Solver(Meshless_Sweep const &wrs);

        // Solve problem
        virtual void solve(std::vector<double> &x) const override;

    protected:
        
        // Data
        std::shared_ptr<Epetra_Comm> comm_;
        std::shared_ptr<Epetra_Map> map_;
        std::vector<std::shared_ptr<Epetra_CrsMatrix> > mat_;
        mutable std::shared_ptr<Epetra_Vector> lhs_;
        mutable std::shared_ptr<Epetra_Vector> rhs_;
        std::vector<std::shared_ptr<Epetra_LinearProblem> > problem_;
        std::vector<std::shared_ptr<Amesos_BaseSolver> > solver_;
    };

    // Amesos parallel solver
    // Solves ordinates in parallel
    class Amesos_Parallel_Solver : public Trilinos_Solver
    {
    public:
        // Constructor
        Amesos_Parallel_Solver(Meshless_Sweep const &wrs);

        // Solve problem
        virtual void solve(std::vector<double> &x) const override;
        
    protected:
        
        // Data
        std::vector<std::shared_ptr<Epetra_Comm> > comm_;
        std::vector<std::shared_ptr<Epetra_Map> > map_;
        std::vector<std::shared_ptr<Epetra_CrsMatrix> > mat_;
        mutable std::vector<std::shared_ptr<Epetra_Vector> > lhs_;
        mutable std::vector<std::shared_ptr<Epetra_Vector> > rhs_;
        std::vector<std::shared_ptr<Epetra_LinearProblem> > problem_;
        std::vector<std::shared_ptr<Amesos_BaseSolver> > solver_;
    };
    
    // Aztec solver
    // Iterative, does not store matrices
    // Only in serial
    class Aztec_Solver : public Trilinos_Solver
    {
    public:
        // Solver options
        struct Options
        {
            int max_iterations = 1000;
            int kspace = 20;
            double tolerance = 1e-8;
        };
        
        // Constructor
        Aztec_Solver(Meshless_Sweep const &wrs);
        
        // Solve problem
        virtual void solve(std::vector<double> &x) const override;

    protected:

        // Data
        std::shared_ptr<Epetra_Comm> comm_;
        std::shared_ptr<Epetra_Map> map_;
        mutable std::shared_ptr<Epetra_Vector> lhs_;
        mutable std::shared_ptr<Epetra_Vector> rhs_;
        
    };

    // Aztec preconditioned by Ifpack
    // Preconditioned by inverse of Linv matrices
    // Only in serial
    class Aztec_Ifpack_Solver : public Trilinos_Solver
    {
    public:
        
        // Constructor
        Aztec_Ifpack_Solver(Meshless_Sweep const &wrs);
        
        // Solve problem
        virtual void solve(std::vector<double> &x) const override;

    protected:
        
        std::shared_ptr<Epetra_Comm> comm_;
        std::shared_ptr<Epetra_Map> map_;
        std::vector<std::shared_ptr<Epetra_CrsMatrix> > mat_;
        mutable std::shared_ptr<Epetra_Vector> lhs_;
        mutable std::shared_ptr<Epetra_Vector> rhs_;
        std::vector<std::shared_ptr<Epetra_LinearProblem> > problem_;
        std::vector<std::shared_ptr<Ifpack_Preconditioner> > prec_;
        std::vector<std::shared_ptr<AztecOO> > solver_;
    };

    // Belos solver: iterative, not preconditioned
    // Works in parallel
    class Belos_Solver : public Trilinos_Solver
    {
    public:
        
        // Constructor
        Belos_Solver(Meshless_Sweep const &wrs);
        
        // Solve problem
        virtual void solve(std::vector<double> &x) const override;
        
    protected:

        std::vector<std::shared_ptr<Epetra_Comm> > comm_;
        std::vector<std::shared_ptr<Epetra_Map> > map_;
        mutable std::vector<std::shared_ptr<Epetra_Vector> > lhs_;
        mutable std::vector<std::shared_ptr<Epetra_Vector> > rhs_;
        std::vector<std::shared_ptr<BelosLinearProblem> > problem_;
        std::vector<std::shared_ptr<BelosSolver> > solver_;
    };

    // Belos preconditioned by Ifpack
    // Preconditioned by inverse of Linv matrices
    // Works in parallel
    class Belos_Ifpack_Solver : public Trilinos_Solver
    {
    public:
        
        // Constructor
        Belos_Ifpack_Solver(Meshless_Sweep const &wrs);
        
        // Solve problem
        virtual void solve(std::vector<double> &x) const override;

    protected:
        
        std::vector<std::shared_ptr<Epetra_Comm> > comm_;
        std::vector<std::shared_ptr<Epetra_Map> > map_;
        std::vector<std::shared_ptr<Epetra_CrsMatrix> > mat_;
        mutable std::vector<std::shared_ptr<Epetra_Vector> > lhs_;
        mutable std::vector<std::shared_ptr<Epetra_Vector> > rhs_;
        std::vector<std::shared_ptr<BelosPreconditioner> > prec_;
        std::vector<std::shared_ptr<BelosLinearProblem> > problem_;
        std::vector<std::shared_ptr<BelosSolver> > solver_;
    };

    // Belos preconditioned on the right by Ifpack
    // Preconditioned on the right by the inverse of the basis value matrix
    // Does not store the matrices between iterations
    // Works in parallel
    class Belos_Ifpack_Right_Solver : public Trilinos_Solver
    {
    public:
        
        // Constructor
        Belos_Ifpack_Right_Solver(Meshless_Sweep const &wrs);
        
        // Solve problem
        virtual void solve(std::vector<double> &x) const override;

    protected:
        
        std::vector<std::shared_ptr<Epetra_Comm> > comm_;
        std::vector<std::shared_ptr<Epetra_Map> > map_;
        std::vector<std::shared_ptr<Epetra_CrsMatrix> > mat_;
        mutable std::vector<std::shared_ptr<Epetra_Vector> > lhs_;
        mutable std::vector<std::shared_ptr<Epetra_Vector> > rhs_;
        std::vector<std::shared_ptr<Epetra_CrsMatrix> > prec_mat_;
        std::vector<std::shared_ptr<BelosPreconditioner> > prec_;
        std::vector<std::shared_ptr<BelosLinearProblem> > problem_;
        std::vector<std::shared_ptr<BelosSolver> > solver_;
    };

    // Belos preconditioned on the right by Ifpack
    // Preconditioned on the right by the inverse of the basis value matrix
    // Stores the matrices between iterations
    // Works in parallel
    class Belos_Ifpack_Right2_Solver : public Trilinos_Solver
    {
    public:
        
        // Constructor
        Belos_Ifpack_Right2_Solver(Meshless_Sweep const &wrs);
        
        // Solve problem
        virtual void solve(std::vector<double> &x) const override;

    protected:
        
        std::vector<std::shared_ptr<Epetra_Comm> > comm_;
        std::vector<std::shared_ptr<Epetra_Map> > map_;
        std::vector<std::shared_ptr<Epetra_CrsMatrix> > mat_;
        mutable std::vector<std::shared_ptr<Epetra_Vector> > lhs_;
        mutable std::vector<std::shared_ptr<Epetra_Vector> > rhs_;
        std::vector<std::shared_ptr<Epetra_CrsMatrix> > prec_mat_;
        std::vector<std::shared_ptr<BelosPreconditioner> > prec_;
        std::vector<std::shared_ptr<BelosLinearProblem> > problem_;
        std::vector<std::shared_ptr<BelosSolver> > solver_;
    };
    
    // Data
    Options options_;
    std::shared_ptr<Weak_Spatial_Discretization> spatial_discretization_;
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;
    std::shared_ptr<Sweep_Solver> solver_;
};

#endif
