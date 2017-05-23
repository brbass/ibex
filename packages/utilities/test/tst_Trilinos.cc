#include <iomanip>
#include <iostream>
#include <string>

#include <mpi.h>

#include <Amesos.h>
#include <AztecOO.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>
#include <Ifpack.h>

#include "Random_Number_Generator.hh"
#include "Timer.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;

// Get random number generator for problem
Random_Number_Generator<double> rng(0, // lower bound
                                    1, // upper bound
                                    581); // seed

// Get all matrices in input file
vector<shared_ptr<Epetra_CrsMatrix> > get_matrices(string filename)
{
    Timer timer;
    timer.start();
    
    // Get XML input document
    XML_Document input_document(filename);
    XML_Node input_node = input_document.get_child("output");

    // Get Epetra comm
    shared_ptr<Epetra_MpiComm> comm
        = make_shared<Epetra_MpiComm>(MPI_COMM_WORLD);

    // Get matrices
    vector<shared_ptr<Epetra_CrsMatrix> > matrices;
    for (XML_Node matrix_node = input_node.get_child("matrix");
         matrix_node;
         matrix_node = matrix_node.get_sibling("matrix",
                                               false))
    {
        // Get size data
        int number_of_points
            = matrix_node.get_attribute<int>("number_of_points");
        vector<int> number_of_entries
            = matrix_node.get_child_vector<int>("number_of_entries",
                                                number_of_points);

        // Get Epetra map
        shared_ptr<Epetra_Map> map
            = make_shared<Epetra_Map>(number_of_points,
                                      0,
                                      *comm);

        // Initialize Epetra matrix
        shared_ptr<Epetra_CrsMatrix> matrix
            = make_shared<Epetra_CrsMatrix>(Copy, // Data access
                                            *map,
                                            &number_of_entries[0],
                                            true); // Static proflile

        // Insert matrix values
        for (XML_Node row_node = matrix_node.get_child("row");
             row_node;
             row_node = row_node.get_sibling("row",
                                             false))
        {
            int index = row_node.get_attribute<int>("row_index");
            vector<int> indices
                = row_node.get_child_vector<int>("column_indices",
                                                 number_of_entries[index]);
            vector<double> values
                = row_node.get_child_vector<double>("values",
                                                    number_of_entries[index]);
            matrix->InsertGlobalValues(index,
                                       number_of_entries[index],
                                       &values[0],
                                       &indices[0]);
        }
        matrix->FillComplete();
        matrix->OptimizeStorage();
        matrices.push_back(matrix);
    }

    timer.stop();
    
    cout << "time to read in matrices: " << timer.time() << endl;
    
    return matrices;
}

shared_ptr<Epetra_Vector> get_random_vector(shared_ptr<Epetra_CrsMatrix> mat)
{
    int number_of_points = mat->NumGlobalRows();
    shared_ptr<Epetra_Vector> vec
        = make_shared<Epetra_Vector>(mat->RowMap());
    vector<double> values = rng.vector(number_of_points);
    for (int i = 0; i < number_of_points; ++i)
    {
        (*vec)[i] = values[i];
    }

    return vec;
}

double get_checksum(shared_ptr<Epetra_Vector> vec)
{
    int number_of_points = vec->GlobalLength();

    double checksum = 0;
    for (int i = 0; i < number_of_points; ++i)
    {
        checksum += (*vec)[i];
    }

    return checksum;
}

void test_lu_decomposition(shared_ptr<Epetra_CrsMatrix> mat,
                           shared_ptr<Epetra_Vector> rhs,
                           double &checksum,
                           double &setup_time,
                           double &solve_time)
{
    // Get data
    int number_of_points = mat->NumGlobalRows();
    shared_ptr<Epetra_Vector> lhs
        = make_shared<Epetra_Vector>(mat->RowMap());
    shared_ptr<Epetra_LinearProblem> problem
        = make_shared<Epetra_LinearProblem>(mat.get(),
                                            lhs.get(),
                                            rhs.get());

    // Setup problem
    Timer timer;
    timer.start();
    Amesos factory;
    shared_ptr<Amesos_BaseSolver*> solver
        = make_shared<Amesos_BaseSolver*>(factory.Create("Klu",
                                                         *problem));
    (*solver)->SymbolicFactorization();
    (*solver)->NumericFactorization();
    timer.stop();
    setup_time = timer.time();
    
    // Solve problem
    timer.start();
    int num_solves = 10;
    for (int i = 0; i < num_solves; ++i)
    {
        (*lhs)[0] = i;
        (*solver)->Solve();
    }
    timer.stop();
    solve_time = timer.time() / num_solves;

    // Get checksum
    checksum = get_checksum(lhs);
}

void test_ilut(shared_ptr<Epetra_CrsMatrix> mat,
                shared_ptr<Epetra_Vector> rhs,
                double &checksum,
                double &setup_time,
                double &solve_time)
{
    // Get data
    int number_of_points = mat->NumGlobalRows();
    shared_ptr<Epetra_Vector> lhs
        = make_shared<Epetra_Vector>(mat->RowMap());
    shared_ptr<Epetra_LinearProblem> problem
        = make_shared<Epetra_LinearProblem>(mat.get(),
                                            lhs.get(),
                                            rhs.get());

    // Setup problem
    Timer timer;
    timer.start();
    int kspace = 20;
    shared_ptr<AztecOO> solver
        = make_shared<AztecOO>(*problem);
    solver->SetAztecOption(AZ_solver, AZ_gmres);
    solver->SetAztecOption(AZ_kspace, kspace);
    solver->SetAztecOption(AZ_precond, AZ_dom_decomp);
    solver->SetAztecOption(AZ_subdomain_solve, AZ_ilut);
    solver->SetAztecOption(AZ_output, AZ_warnings);
    solver->SetAztecOption(AZ_keep_info, 1); // Keeps info for multiple solves
    double condest;
    solver->ConstructPreconditioner(condest);
    if (condest > 1e14)
    {
        cout << "preconditioner condition number too large" << endl;
    }
    solver->SetAztecOption(AZ_pre_calc, AZ_reuse); // Prevents recomputation of preconditioner
    timer.stop();
    setup_time = timer.time();
    
    // Solve problem
    int max_iterations = 100;
    double tolerance = 1e-8;
    int num_solves = 10;
    timer.start();
    for (int i = 0; i < num_solves; ++i)
    {
        (*lhs)[0] = i;
        solver->Iterate(max_iterations,
                        tolerance);
    }
    timer.stop();
    solve_time = timer.time() / num_solves;
    
    // Get checksum
    checksum = get_checksum(lhs);
}

void test_ifpack(shared_ptr<Epetra_CrsMatrix> mat,
                 shared_ptr<Epetra_Vector> rhs,
                 double &checksum,
                 double &setup_time,
                 double &solve_time)
{
    // Get data
    int number_of_points = mat->NumGlobalRows();
    shared_ptr<Epetra_Vector> lhs
        = make_shared<Epetra_Vector>(mat->RowMap());
    shared_ptr<Epetra_LinearProblem> problem
        = make_shared<Epetra_LinearProblem>(mat.get(),
                                            lhs.get(),
                                            rhs.get());

    // Setup problem
    Timer timer;
    timer.start();
    Ifpack factory;
    shared_ptr<Ifpack_Preconditioner*> prec
        = make_shared<Ifpack_Preconditioner*>(factory.Create("ILU",
                                                             mat.get()));
    Teuchos::ParameterList list;
    // list.set("fact: drop tolerance", 1e-10);
    list.set("fact: level-of-fill", 2);
    (*prec)->SetParameters(list);
    (*prec)->Initialize();
    (*prec)->Compute();
    int kspace = 10;
    shared_ptr<AztecOO> solver
        = make_shared<AztecOO>(*problem);
    solver->SetAztecOption(AZ_solver, AZ_gmres);
    solver->SetAztecOption(AZ_kspace, kspace);
    solver->SetPrecOperator(*prec);
    solver->SetAztecOption(AZ_output, AZ_warnings);
    timer.stop();
    setup_time = timer.time();
    
    // Solve problem
    int max_iterations = 100;
    double tolerance = 1e-8;
    int num_solves = 10;
    timer.start();
    for (int i = 0; i < num_solves; ++i)
    {
        (*lhs)[0] = i;
        solver->Iterate(max_iterations,
                        tolerance);
    }
    timer.stop();
    solve_time = timer.time() / num_solves;
    
    // Get checksum
    checksum = get_checksum(lhs);
}

void run_test(shared_ptr<Epetra_CrsMatrix> mat)
{
    // Get RHS vector
    int number_of_points = mat->NumGlobalRows();
    shared_ptr<Epetra_Vector> rhs
        = get_random_vector(mat);

    // Get lists for methods
    int num_methods = 3;
    vector<string> descriptions
        = {"full lu",
           "ilut",
           "ifpack"};
    vector<function<void(shared_ptr<Epetra_CrsMatrix>,
                         shared_ptr<Epetra_Vector>,
                         double &,
                         double &,
                         double &)> > methods
        = {test_lu_decomposition,
           test_ilut,
           test_ifpack};

    // Run methods
    int w = 16;
    cout << setw(w) << "method";
    cout << setw(w) << "checksum";
    cout << setw(w) << "setup_time";
    cout << setw(w) << "solve_time";
    cout << endl;
    for (int i = 0; i < num_methods; ++i)
    {
        double checksum;
        double setup_time;
        double solve_time;
        methods[i](mat,
                   rhs,
                   checksum,
                   setup_time,
                   solve_time);
        cout << setw(w) << descriptions[i];
        cout << setw(w) << checksum;
        cout << setw(w) << setup_time;
        cout << setw(w) << solve_time;
        cout << endl;
    }
    cout << endl;
}

void run_tests(string filename)
{
    cout << "running tests" << endl;
    vector<shared_ptr<Epetra_CrsMatrix> > matrices
        = get_matrices(filename);
    for (shared_ptr<Epetra_CrsMatrix> matrix : matrices)
    {
        run_test(matrix);
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    if (argc != 2)
    {
        cerr << "usage: tst_trilinos [input.xml]" << endl;
        return 1;
    }
    
    string filename = argv[1];
    
    run_tests(filename);

    MPI_Finalize();
}
