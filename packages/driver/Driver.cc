#include "Driver.hh"

#if defined(ENABLE_OPENMP)
    #include <omp.h>
#else
    inline void omp_set_num_threads(int i) {return;}
#endif

#include "Integration_Problem.hh"
#include "Heat_Problem.hh"
#include "Manufactured_Problem.hh"
#include "Transport_Problem.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using std::string;

Driver::
Driver(string input_filename):
    input_filename_(input_filename)
{
    output_filename_ = input_filename + ".out";
    
    run_problem();    
}

void Driver::
run_problem()
{
    // Get input file
    XML_Document input_file(input_filename_);
    XML_Node input_node = input_file.get_child("input");

    // Get output file
    XML_Document output_file;
    XML_Node output_node = output_file.append_child("output");

    // Get problem type
    string problem_type = input_node.get_attribute<string>("type");
    bool print_progress = input_node.get_attribute<bool>("print",
                                                         false);
    int number_of_threads = input_node.get_attribute<int>("number_of_threads",
                                                    1);
    omp_set_num_threads(number_of_threads);
    
    // Solve problem
    if (problem_type == "transport")
    {
        Transport_Problem transport_problem(input_node,
                                            output_node,
                                            print_progress);
        transport_problem.solve();
    }
    else if (problem_type == "manufactured")
    {
        Manufactured_Problem manufactured_problem(input_node,
                                                  output_node,
                                                  print_progress);
        manufactured_problem.solve();
    }
    else if (problem_type == "heat")
    {
        Heat_Problem heat_problem(input_node,
                                  output_node,
                                  print_progress);
        heat_problem.solve();
    }
    else if (problem_type == "integration")
    {
        Integration_Problem int_problem(input_node,
                                        output_node,
                                        print_progress);
        int_problem.solve();
    }
    else
    {
        AssertMsg(false, "problem type (" + problem_type + ") not found");
    }

    // Save output file
    output_file.save(output_filename_);
}

