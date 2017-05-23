#ifndef Weak_RBF_Sweep_hh
#define Weak_RBF_Sweep_hh

#include "Sweep_Operator.hh"
#include "Weak_Spatial_Discretization.hh"

class Amesos_BaseSolver;
class AztecOO;
class Epetra_CrsMatrix;
class Epetra_Comm;
class Epetra_LinearProblem;
class Epetra_Map;
class Epetra_Vector;
class Ifpack_Preconditioner;

class Weak_RBF_Sweep : public Sweep_Operator
{
public:

    // Sweep options
    struct Options
    {
        enum class Solver
        {
            AMESOS,
            AZTEC,
            AZTEC_ILUT,
            AZTEC_IFPACK
        };
        
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
    virtual void check_class_invariants() const override;
    virtual std::string description() const override
    {
        return "Weak_RBF_Sweep";
    }

    // Save matrix to specified XML output file
    void save_matrix_as_xml(int o,
                            int g,
                            XML_Node output_node) const;
    
private:

    // Vector_Operator function
    virtual void apply(std::vector<double> &x) const override;

    // Weak_RBF_Sweep functions
    void update_augments(std::vector<double> &x) const;
    void get_matrix_row(int i, // weight function index (row)
                        int o, // ordinate
                        int g, // group
                        std::vector<int> &indices, // global basis (column indices)
                        std::vector<double> &values) const; // column values
    void get_rhs(int i, // weight function index (row)
                 int o, // ordinate
                 int g, // group
                 std::vector<double> const &x, // angular flux w/ augments
                 double &value) const; // rhs value
    
    // Generalized solver
    class Sweep_Solver
    {
    public:
        // Constructor
        Sweep_Solver(Weak_RBF_Sweep const &wrs);

        // Solve problem
        virtual void solve(std::vector<double> &x) const = 0;

    protected:
        // Data
        Weak_RBF_Sweep const &wrs_;
    };
    
    // Generalized trilinos solver
    class Trilinos_Solver : public Sweep_Solver
    {
    public:
        // Constructor
        Trilinos_Solver(Weak_RBF_Sweep const &wrs);

        // Solve problem
        virtual void solve(std::vector<double> &x) const = 0;

    protected:
        
        // Return transport matrix for o and g
        std::shared_ptr<Epetra_CrsMatrix> get_matrix(int o,
                                                     int g) const;

        // Set rhs_ to current source
        void set_rhs(int o,
                     int g,
                     std::vector<double> const &x) const;
        
        // Data
        std::shared_ptr<Epetra_Comm> comm_;
        std::shared_ptr<Epetra_Map> map_;
        mutable std::shared_ptr<Epetra_Vector> lhs_;
        mutable std::shared_ptr<Epetra_Vector> rhs_;
    };
    
    // Amesos solver: stores LU decompositions of all matrices
    class Amesos_Solver : public Trilinos_Solver
    {
    public:
        // Constructor
        Amesos_Solver(Weak_RBF_Sweep const &wrs);

        // Solve problem
        virtual void solve(std::vector<double> &x) const override;

    protected:
        // Data
        std::vector<std::shared_ptr<Epetra_CrsMatrix> > mat_;
        std::vector<std::shared_ptr<Epetra_LinearProblem> > problem_;
        std::vector<std::shared_ptr<Amesos_BaseSolver> > solver_;
    };
    
    // Aztec solver: iterative, does not store matrices
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
        Aztec_Solver(Weak_RBF_Sweep const &wrs,
                     Options options);
        
        // Solve problem
        virtual void solve(std::vector<double> &x) const override;

    protected:
        Options options_;
    };
    
    // Aztec solver: iterative, stores partial LU decompositions
    // Not working: does not find past preconditioner
    class Aztec_ILUT_Solver : public Trilinos_Solver
    {
    public:
        // Solver options
        struct Options
        {
            int max_iterations = 100;
            int kspace = 20;
            double tolerance = 1e-8;
        };
        
        // Constructor
        Aztec_ILUT_Solver(Weak_RBF_Sweep const &wrs,
                          Options options);
        
        // Solve problem
        virtual void solve(std::vector<double> &x) const override;

    protected:
        Options options_;
        std::vector<std::shared_ptr<Epetra_CrsMatrix> > mat_;
        std::vector<std::shared_ptr<Epetra_LinearProblem> > problem_;
        std::vector<std::shared_ptr<AztecOO> > solver_;
    };

    class Aztec_Ifpack_Solver : public Trilinos_Solver
    {
    public:
        // Solver options
        struct Options
        {
            int max_iterations = 100;
            int kspace = 20;
            int level_of_fill = 2;
            double tolerance = 1e-8;
            // double drop_tolerance = 1e-10;
        };
        
        // Constructor
        Aztec_Ifpack_Solver(Weak_RBF_Sweep const &wrs,
                            Options options);
        
        // Solve problem
        virtual void solve(std::vector<double> &x) const override;

    protected:
        Options options_;
        std::vector<std::shared_ptr<Epetra_CrsMatrix> > mat_;
        std::vector<std::shared_ptr<Epetra_LinearProblem> > problem_;
        std::vector<std::shared_ptr<Ifpack_Preconditioner> > prec_;
        std::vector<std::shared_ptr<AztecOO> > solver_;
    };
    
    // Data
    Options options_;
    std::shared_ptr<Weak_Spatial_Discretization> spatial_discretization_;
    std::shared_ptr<Angular_Discretization> angular_discretization_;
    std::shared_ptr<Energy_Discretization> energy_discretization_;
    std::shared_ptr<Sweep_Solver> solver_;
};

#endif
