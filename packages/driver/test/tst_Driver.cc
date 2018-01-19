#include <iostream>
#include <string>

#include <mpi.h>

#include "Check_Equality.hh"
#include "Driver.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;

int main(int argc, char **argv)
{
    int checksum = 0;
    
    MPI_Init(&argc, &argv);
    if (argc != 2)
    {
        cerr << "usage: tst_driver [input.xml]" << endl;
        return 1;
    }

    // Get input, output and result files
    string input_filename = argv[1];
    Driver driver(input_filename);
    string output_filename = input_filename + ".out";
    string result_filename = input_filename + ".res";

    XML_Document output_document(output_filename);
    XML_Document result_document(result_filename);

    XML_Node output_node = output_document.get_child("output");
    XML_Node result_cases = result_document.get_child("results");
    
    // Compare results for each given case
    for (XML_Node result_node = result_cases.get_child("result");
         result_node;
         result_node = result_node.get_sibling("result",
                                               false))
    {
        int expected_size = result_node.get_attribute<int>("size");
        double tolerance = result_node.get_attribute<double>("tolerance");
        vector<string> comparison_location = result_node.get_attribute_vector<string>("location");
        XML_Node output_child = output_node;
        for (string location : comparison_location)
        {
            output_child = output_child.get_child(location);
        }
        vector<double> result_data = result_node.get_child_vector<double>("data",
                                                                          expected_size);
        vector<double> output_data = output_child.get_vector<double>(expected_size);
        
        if (!Check_Equality::approx(result_data, output_data, tolerance))
        {
            cout << "result in node (" << output_child.name() << ") incorrect" << endl;
            checksum += 1;
        }
    }
    
    MPI_Finalize();
    
    return checksum;
}

