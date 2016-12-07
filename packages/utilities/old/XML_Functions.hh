#ifndef XML_Functions_hh
#define XML_Functions_hh

#include <algorithm>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "pugixml.hh"

#include "Check.hh"
#include "Vector_String_Functions.hh"

namespace XML_Functions
{
    using std::setprecision;
    using std::string;
    using std::vector;
    
    /*
      Functions to make using the XML parser more intuitive
    */
    
    // Get a child value of the node
    template<typename T> T 
    child_value(pugi::xml_node &node, 
                string description,
                bool required = true);
    
    // Get a vector of values from the node
    template<typename T> vector<T> 
    child_vector(pugi::xml_node &node,
                 string description, 
                 unsigned expected_size,
                 bool required = true)
    {
        vector<T> value;
    
        string_to_vector(child_value<string>(node, description, required),
                         value);
    
        AssertMsg(value.size() == expected_size, description + " size");
    
        return value;
    }
    
    // Append a child to the XML file
    template<typename T> void
    append_child(pugi::xml_node &node,
                 T data,
                 string description,
                 unsigned precision = 0)
    {
        std::stringstream ss;
        ss << setprecision(XML_PRECISION);
        ss << data;
        
        string data_string = ss.str();
        
        pugi::xml_node child = node.append_child(description.c_str());
        
        child.append_child(pugi::node_pcdata).set_value(data_string.c_str());
    }

    // Get attribute
    template<typename T> T
    attribute_value(pugi::xml_node &node,
                    string description,
                    bool required = true);
    
    // Specializations
    template<> inline string
    child_value<string> (pugi::xml_node &node,
                         string description,
                         bool required)
    {
        string value = static_cast<string>(node.child_value(description.c_str()));

        if (required && value == "")
        {
            string name = static_cast<string>(node.name());
            string error_message
                = "required value ("
                + description
                + ") in node ("
                + name
                + ") not found";
            
            AssertMsg(false, error_message);
        }
    
        return value;
    }

    template<> inline int
    child_value<int> (pugi::xml_node &node, 
                      string description,
                      bool required)
    {
        return stoi(child_value<string>(node, description, required));
    }

    template<> inline double
    child_value<double> (pugi::xml_node &node,
                         string description,
                         bool required)
    {
        return stof(child_value<string>(node, description, required));
    }
    
    template<typename T> void
    append_child(pugi::xml_node &node,
                 vector<T> const &data,
                 string description,
                 string index_order = "")
    {
        string data_string;
        vector_to_string(data_string,
                         data);
    
        pugi::xml_node child = node.append_child(description.c_str());
    
        child.append_child(pugi::node_pcdata).set_value(data_string.c_str());

        if (index_order != "")
        {
            pugi::xml_attribute it = child.append_attribute("data_order");
            it.set_value(index_order.c_str());
        }
    }
    
    // template<> int
    // attribute_value(pugi::xml_node &node,
    //                 string description,
    //                 bool required)
    // {
    //     pugi::xml_attribute attr = node.attribute(description.c_str());
    //     if (required && attr.empty())
    //     {
    //         string name = static_cast<string>(node.name());
    //         string error_message
    //             = "required value ("
    //             + description
    //             + ") in node ("
    //             + name
    //             + ") not found";
            
    //         AssertMsg(false, error_message);
    //     }
        
    //     return attr.as_int()
    // }
} // namespace XML_Functions

#endif
