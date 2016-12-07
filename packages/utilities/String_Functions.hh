#ifndef String_Functions_hh
#define String_Functions_hh

#include <algorithm>
#include <iomanip>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace String_Functions
{
    using std::string;
    using std::vector;

    // Convert data to string
    template<class T> string
    to_string(string &data_string,
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
    string_to_vector(string const &data_string,
                     vector<T> &data)
    {
        std::istringstream iss(data_string);
        
        data = vector<T>{std::istream_iterator<T>(iss), std::istream_iterator<T>()};
    }

    // Convert a vector to a string
    template<class T> void
    vector_to_string(string &data_string,
                     vector<T> const &data,
                     int precision = 16)
    {
        std::ostringstream oss;
        oss << std::setprecision(precision);
        std::copy(data.begin(), data.end(), std::ostream_iterator<T>(oss, " "));
        
        data_string = oss.str();
    }
}

#endif
