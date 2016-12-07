#ifndef XML_Functions_New_hh
#define XML_Functions_New_hh

#include <string>
#include <vector>

#include "pugixml.hh"

#define XML_PRECISION 16

namespace XML_Functions
{
    // Convert a string to a vector
    template<class T> void 
    string_to_vector(string const &data_string,
                     vector<T> &data);

    // Convert a vector to a string
    template<class T> void
    vector_to_string(string &data_string,
                     vector<T> const &data);
    
    // Function definitions
    template<class T> void 
    string_to_vector(string const &data_string,
                     vector<T> &data)
    {
        std::istringstream iss(data_string);
        
        data = vector<T>{std::istream_iterator<T>(iss), std::istream_iterator<T>()};
    }
    
    template<class T> void
    vector_to_string(string &data_string,
                     vector<T> const &data)
    {
        std::ostringstream oss;
        oss << setprecision(XML_PRECISION);
        std::copy(data.begin(), data.end(), std::ostream_iterator<T>(oss, " "));
        
        data_string = oss.str();
    }
    
    // Get child value
    template<typename T> T child_value(pugi::xml_node &value_node);
    
    // Child value with default value
    template<typename T> T child_value_def(pugi::xml_node &node,
                                           string desc,
                                           T &def)
    {
        pugi::xml_node child = node.child(desc);

        if (child.empty())
        {
            return def;
        }
    }

    
}

#endif
