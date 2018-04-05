#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "Quadrature_Rule.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;

void output_quadrature(XML_Node input_node,
                       XML_Node output_node)
{
    string geometry = input_node.get_child_value<string>("geometry");

    Quadrature_Rule::Quadrature_Type quad_type = Quadrature_Rule::Quadrature_Type::GAUSS_LEGENDRE;
    if (geometry == "double_cylindrical")
    {
        vector<int> num = input_node.get_child_vector<int>("num", 2);
        vector<double> centers = input_node.get_child_vector<double>("centers", 4);
        vector<double> radii = input_node.get_child_vector<double>("radii", 2);
        vector<double> boundaries = input_node.get_child_vector<double>("boundaries", 4);
        vector<double> ordinates_x;
        vector<double> ordinates_y;
        vector<double> weights;
        Quadrature_Rule::cartesian_bounded_double_cylindrical_2d(quad_type,
                                                                 quad_type,
                                                                 num[0],
                                                                 num[1],
                                                                 centers[0],
                                                                 centers[1],
                                                                 radii[0],
                                                                 centers[2],
                                                                 centers[3],
                                                                 radii[1],
                                                                 boundaries[0],
                                                                 boundaries[1],
                                                                 boundaries[2],
                                                                 boundaries[3],
                                                                 ordinates_x,
                                                                 ordinates_y,
                                                                 weights);
        int num_ordinates = weights.size();
        output_node.set_child_value(num_ordinates, "number_of_ordinates");
        output_node.set_child_vector(ordinates_x, "ordinates_x");
        output_node.set_child_vector(ordinates_y, "ordinates_y");
        output_node.set_child_vector(weights, "weights");
    }
    else
    {
        cout << "quad type not implemented" << endl;
    }
}

int main(int argc, char **argv)
{
    string input_filename = argv[1];
    XML_Document input_file(input_filename);
    XML_Node input_node = input_file.get_child("input");
    string output_filename = input_filename + ".out";
    XML_Document output_file;
    XML_Node output_node = output_file.append_child("output");
    output_quadrature(input_node,
                      output_node);
    output_file.save(output_filename);
    
    return 0;
}
