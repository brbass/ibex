#include <iostream>
#include <string>

#include "Check_Equality.hh"
#include "XML_Document.hh"
#include "XML_Node.hh"

using namespace std;

namespace ce = Check_Equality;

int test_xml_read(string input_folder)
{
    int checksum = 0;
    
    string input_filename = input_folder + "/test_xml.xml";

    XML_Document doc(input_filename);
    XML_Node input = doc.get_child("input");
    XML_Node cat = input.get_child("cat");

    // Test attributes
    string breed = cat.get_attribute_value<string>("breed");
    if (breed != "siamese")
    {
        cerr << "tst_XML: breed incorrect" << endl;
        checksum += 1;
    }

    // Test scalars
    string color = cat.get_child_value<string>("color");
    if (color != "black")
    {
        cerr << "tst_XML: color incorrect" << endl;
        checksum += 1;
    }
    XML_Node bites = cat.get_child("bites");
    if (!bites.get_value<bool>())
    {
        cerr << "tst_XML: all cats bite" << endl;
    }

    // Test vectors
    vector<double> dimensions = cat.get_child_vector<double>("dimensions",
                                                             3);
    if (!ce::approx(dimensions, {14.2, 5.9, 9.0}, 1e-14))
    {
        cerr << "tst_XML: dimensions incorrect" << endl;
        checksum += 1;
    }

    // Test default values
    int hooves = cat.get_attribute_value<int>("hooves",
                                              0);
    if (hooves != 0)
    {
        cerr << "tst_XML: cats do not have hooves" << endl;
        checksum += 1;
    }
    bool is_alive = cat.get_value<bool>(true);
    if (!is_alive)
    {
        cerr << "tst_XML: soulless != dead" << endl;
        checksum += 1;
    }
    vector<int> lives_left = cat.get_vector<int>(3, {1, 2, 4});
    if (!ce::equal(lives_left, {1, 2, 4}))
    {
        cerr << "tst_XML: cat has died 6 times" << endl;
        checksum += 1;
    }
    
    return checksum;
}

int main(int argc, char **argv)
{
    int checksum = 0;

    if (argc != 2)
    {
        cerr << "usage: tst_XML [input_folder]" << endl;
        return 1;
    }
    
    string input_folder = argv[1];
    
    checksum += test_xml_read(input_folder);
    
    return checksum;
}
