#ifndef String_Functions_hh
#define String_Functions_hh

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

namespace String_Functions
{
    // Convert data to string
    template<class T> std::string
    to_string(std::string &data_string,
              T &data,
              int precision = 16)
    {
        std::stringstream ss;
        ss << std::setprecision(precision);
        ss << data;
        
        return ss.str();
    }
    
    // Convert a string to a vector
    template<class T> void 
    string_to_vector(std::string const &data_string,
                     std::vector<T> &data)
    {
        std::istringstream iss(data_string);
        
        data = std::vector<T>{std::istream_iterator<T>(iss), std::istream_iterator<T>()};
    }

    // Convert a vector to a string
    template<class T> void
    vector_to_string(std::string &data_string,
                     std::vector<T> const &data,
                     int precision = 16)
    {
        std::ostringstream oss;
        oss << std::setprecision(precision);
        std::copy(data.begin(), data.end(), std::ostream_iterator<T>(oss, " "));
        
        data_string = oss.str();
    }
}

#endif
